#!/usr/bin/env python3

import argparse
import csv
import json
import fcntl
import os
import re
import socket
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime

import pandas as pd


NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_URL = "https://rest.uniprot.org/taxonomy/{taxid}"

EXCLUDED_RANKS = {"no rank", "clade"}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate taxonomy metadata.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input TSV.")
    parser.add_argument("-o", "--output", default=Path("output"), type=Path, help="Output directory.")
    parser.add_argument("--taxonomy-cache", type=Path, help="Reuse validated cached lineages and mnemonics without remote requests.")
    parser.add_argument("--response-cache", type=Path, help="Persistent cache of validated taxonomy and mnemonic responses, shared across runs.")
    return parser.parse_args()


def retry_after_seconds(value):
    """Interpret both HTTP Retry-After formats without shortening server delays."""
    try:
        return max(0.0, float(value))
    except (ValueError, TypeError):
        try:
            date = parsedate_to_datetime(value)
            return max(0.0, (date - datetime.now(timezone.utc)).total_seconds())
        except (ValueError, TypeError, OverflowError):
            return 0.0


_last_ncbi_request = None


def throttle_ncbi(directory=None):
    """Space requests by one second, also across runs sharing the same cache."""
    global _last_ncbi_request
    if directory is None:
        now = time.monotonic()
        if _last_ncbi_request is not None:
            time.sleep(max(0.0, 1.0 - (now - _last_ncbi_request)))
        _last_ncbi_request = time.monotonic()
        return
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / '.ncbi-rate.lock').open('a+') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.seek(0)
        previous = handle.read().strip()
        delay = max(0.0, 1.0 - (time.time() - float(previous))) if previous else 0.0
        time.sleep(delay)
        handle.seek(0)
        handle.truncate()
        handle.write(str(time.time()))
        handle.flush()


def fetch(url, accept, *, attempts=8, timeout=60, rate_directory=None):
    """Retry rate limiting and temporary failures; never retry permanent HTTP errors."""
    request = urllib.request.Request(
        url,
        headers={
            "Accept": accept,
            "User-Agent": "AnnoGlyph",
        },
    )

    host = urllib.parse.urlparse(url).hostname
    for attempt in range(attempts):
        if host == 'eutils.ncbi.nlm.nih.gov':
            throttle_ncbi(rate_directory)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code not in (408, 429, 500, 502, 503, 504) or attempt == attempts - 1:
                raise
            delay = max(min(120, 2 ** (attempt + 1)),
                        retry_after_seconds(error.headers.get('Retry-After')) if error.headers else 0)
            error.close()
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError):
            if attempt == attempts - 1:
                raise
            delay = min(120, 2 ** (attempt + 1))
        print(f'Temporary taxonomy service failure ({host}); retry {attempt + 2}/{attempts} in {delay:g}s', file=sys.stderr, flush=True)
        time.sleep(delay)


def fetch_taxonomy(taxid, rate_directory=None):
    """Return NCBI taxonomy information."""
    query = urllib.parse.urlencode(
        {
            "db": "taxonomy",
            "id": taxid,
            "retmode": "xml",
            "tool": "AnnoGlyph",
        }
    )

    root = ET.fromstring(fetch(f"{NCBI_EFETCH_URL}?{query}", "application/xml", rate_directory=rate_directory))
    taxon = root.find("Taxon")
    if taxon is None:
        raise ValueError(f'NCBI returned no taxonomy record for taxid {taxid}: {root.findtext("ERROR") or root.tag}')
    lineage = []

    for node in taxon.findall("./LineageEx/Taxon") + [taxon]:
        rank = node.findtext("Rank")

        if rank in EXCLUDED_RANKS:
            continue

        lineage.append(
            {
                "rank": rank,
                "taxid": node.findtext("TaxId"),
                "name": node.findtext("ScientificName"),
            }
        )

    return {
        "taxid": taxon.findtext("TaxId"),
        "scientific_name": taxon.findtext("ScientificName"),
        "lineage": lineage,
    }


def cached_taxonomy(directory, genome_id, taxid):
    """Read the historical lineage for this exact sample and terminal taxon."""
    lineage_path = directory / genome_id / f"{taxid}_lineage.csv"
    mnemonic_path = directory / genome_id / f"{taxid}.mnemonic"
    with lineage_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or rows[-1]["id"] != str(taxid):
        raise ValueError(f"Cached lineage does not end at requested taxon {taxid}: {lineage_path}")
    lineage = [{"rank": row["rank"], "taxid": row["id"], "name": row["name"]}
               for row in rows if row["rank"] not in EXCLUDED_RANKS]
    if len({row["rank"] for row in lineage}) != len(lineage):
        raise ValueError(f"Repeated ranked taxa in cached lineage: {lineage_path}")
    taxonomy = {"taxid": str(taxid), "scientific_name": rows[-1]["name"], "lineage": lineage}
    mnemonic = mnemonic_path.read_text().strip() if mnemonic_path.is_file() else "NA"
    return taxonomy, mnemonic or "NA"


def fetch_mnemonic(taxid):
    """Return the UniProt mnemonic."""
    data = json.loads(fetch(UNIPROT_URL.format(taxid=taxid), "application/json"))
    return str(data.get("mnemonic") or "NA")


def cached_remote(directory, kind, taxid, loader):
    """Persist each successful response immediately, even if a later taxon fails."""
    if directory is None:
        return loader()
    if not str(taxid).isdigit():
        raise ValueError(f'Invalid taxid: {taxid!r}')
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f'{kind}_{taxid}.json'
    with target.with_suffix('.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if target.exists():
            record = json.loads(target.read_text())
            if record.get('taxid') != str(taxid) or record.get('kind') != kind:
                raise ValueError(f'Invalid taxonomy cache identity: {target}')
            return record['value']
        value = loader()
        if kind == 'taxonomy' and (not value.get('lineage') or not value.get('scientific_name')):
            raise ValueError(f'Incomplete taxonomy response for {taxid}')
        with tempfile.NamedTemporaryFile(mode='w', dir=directory, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump({'taxid': str(taxid), 'kind': kind, 'value': value}, handle)
        try:
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)
        return value


def fallback_mnemonic(genome_id):
    """Generate a mnemonic from the genome identifier."""
    parts = genome_id.split("-")

    if len(parts) >= 2:
        return (parts[0][:3] + parts[1][:2]).upper()

    return re.sub(r"[^A-Za-z0-9]", "", genome_id)[:5].upper()


def assign_mnemonics(records):
    """Assign official and unique AnnoGlyph mnemonics."""
    bases = []

    for record in records:
        user = record["user_mnemonic"]
        official = record["official_mnemonic"]

        if user and user != "NA":
            base = user
            official = user
        elif official != "NA" and not official.startswith("9"):
            base = official
        else:
            base = fallback_mnemonic(record["id"])

        record["official_mnemonic"] = official
        record["_base"] = base
        bases.append(base)

    counts = Counter(bases)
    indices = Counter()

    for record in records:
        base = record.pop("_base")

        if counts[base] == 1:
            record["annoglyph_mnemonic"] = base
        else:
            indices[base] += 1
            record["annoglyph_mnemonic"] = f"{base}{indices[base]}"


def build_tree(records):
    """Build a taxonomy tree."""
    tree = {}

    for record in records:
        node = tree

        for taxon in record["taxonomy"]["lineage"]:
            node = node.setdefault(
                taxon["taxid"],
                {
                    "name": taxon["name"],
                    "children": {},
                },
            )["children"]

    return tree


def write_newick(nodes):
    """Convert a taxonomy tree to Newick."""
    children = []

    for node in nodes.values():
        subtree = write_newick(node["children"])

        if subtree:
            children.append(f"({subtree}){node['name']}")
        else:
            children.append(node["name"])

    return ",".join(children)


def write_minimal_newick(nodes):
    """Convert a taxonomy tree to Newick and collapse single-child nodes."""
    children = []

    for node in nodes.values():
        node_children = node["children"]

        if not node_children:
            children.append(node["name"])
            continue

        subtree = write_minimal_newick(node_children)

        if len(node_children) == 1:
            children.append(subtree)
        else:
            children.append(f"({subtree}){node['name']}")

    return ",".join(children)


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    dataframe = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")

    taxonomy_cache = {}
    mnemonic_cache = {}
    records = []
    used_ranks = []

    for _, row in dataframe.iterrows():
        taxid = row["taxid"]

        if args.taxonomy_cache and (args.taxonomy_cache / row['id'] / f'{taxid}_lineage.csv').is_file():
            taxonomy_cache[taxid], mnemonic_cache[taxid] = cached_taxonomy(args.taxonomy_cache, row["id"], taxid)
        else:
            if taxid not in taxonomy_cache:
                taxonomy_cache[taxid] = cached_remote(args.response_cache, 'taxonomy', taxid,
                    lambda: fetch_taxonomy(taxid, args.response_cache))
            if taxid not in mnemonic_cache:
                mnemonic_cache[taxid] = cached_remote(args.response_cache, 'mnemonic', taxid,
                    lambda: fetch_mnemonic(taxid))

        taxonomy = taxonomy_cache[taxid]

        for taxon in taxonomy["lineage"]:
            if taxon["rank"] not in used_ranks:
                used_ranks.append(taxon["rank"])

        records.append(
            {
                "id": row["id"],
                "name": row["name"],
                "taxonomy": taxonomy,
                "user_mnemonic": row.get("mnemonic", ""),
                "official_mnemonic": mnemonic_cache[taxid],
            }
        )

    assign_mnemonics(records)

    rows = []
    taxonomy_json = {}

    for record in records:
        taxonomy = record["taxonomy"]

        row = {
            "id": record["id"],
            "name": record["name"],
            "taxonomy": taxonomy["scientific_name"],
            "taxid": taxonomy["taxid"],
            "official_mnemonic": record["official_mnemonic"],
            "annoglyph_mnemonic": record["annoglyph_mnemonic"],
        }

        lineage = {}

        for taxon in taxonomy["lineage"]:
            row[taxon["rank"]] = taxon["name"]
            lineage[taxon["rank"]] = {
                "name": taxon["name"],
                "taxid": taxon["taxid"],
            }

        rows.append(row)

        taxonomy_json[record["id"]] = {
            "id": record["id"],
            "name": record["name"],
            "taxonomy": taxonomy["scientific_name"],
            "taxid": taxonomy["taxid"],
            "mnemonic": {
                "official": record["official_mnemonic"],
                "annoglyph": record["annoglyph_mnemonic"],
            },
            "lineage": lineage,
        }

    columns = [
        "id",
        "name",
        "taxonomy",
        "taxid",
        "official_mnemonic",
        "annoglyph_mnemonic",
        *used_ranks,
    ]

    pd.DataFrame(rows).reindex(columns=columns).to_csv(
        args.output / "taxonomy.tsv",
        sep="\t",
        index=False,
    )

    with open(args.output / "taxonomy.json", "w", encoding="utf-8") as output:
        json.dump(taxonomy_json, output, indent=2, ensure_ascii=False)
        output.write("\n")

    tree = build_tree(records)

    with open(args.output / "tree.nwk", "w", encoding="utf-8") as output:
        output.write(f"({write_newick(tree)});\n")

    with open(args.output / "tree_minimal.nwk", "w", encoding="utf-8") as output:
        output.write(f"({write_minimal_newick(tree)});\n")


if __name__ == "__main__":
    main()
