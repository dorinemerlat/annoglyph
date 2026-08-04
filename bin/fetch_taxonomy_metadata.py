#!/usr/bin/env python3
"""
Fetch taxonomy metadata for an NCBI Taxonomy identifier.

The script retrieves:
- the ranked taxonomic lineage from NCBI Taxonomy;
- the UniProt mnemonic when available.

Outputs:
- <taxid>_lineage.tsv
- <taxid>.mnemonic

Usage:
    fetch_taxonomy_metadata.py <taxid>
"""

import argparse
import csv
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_TAXONOMY_URL = "https://rest.uniprot.org/taxonomy/{taxid}"

EXCLUDED_RANKS = {"clade", "no rank"}
MAX_RETRIES = 5

def parse_args():
    parser = argparse.ArgumentParser(description="Fetch NCBI taxonomy lineage and UniProt mnemonic.")
    parser.add_argument("taxid", help="NCBI Taxonomy identifier.")
    return parser.parse_args()


def open_url(url, accept, timeout=30):
    """Open a URL with automatic retries on HTTP 429."""
    request = urllib.request.Request(
        url,
        headers={
            "Accept": accept,
            "User-Agent": "EXOGAP taxonomy metadata fetcher",
        },
    )

    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()

        except urllib.error.HTTPError as error:
            if error.code == 429 and attempt < MAX_RETRIES - 1:
                wait = 2 ** attempt
                print(
                    f"HTTP 429 received. Waiting {wait} s before retry "
                    f"({attempt + 1}/{MAX_RETRIES})...",
                    file=sys.stderr,
                )
                time.sleep(wait)
                continue

            response_body = error.read().decode("utf-8", errors="replace")
            raise RuntimeError(
                f"HTTP {error.code} while requesting {url}: "
                f"{error.reason}. Response: {response_body}"
            ) from error

        except urllib.error.URLError as error:
            raise RuntimeError(
                f"Failed to request {url}: {error.reason}"
            ) from error

    raise RuntimeError(f"Unable to retrieve {url}")


def parse_ncbi_taxon(taxon_element):
    """Convert an NCBI Taxon XML element into a lineage entry."""
    rank = taxon_element.findtext("Rank", default="").strip()
    taxid = taxon_element.findtext("TaxId", default="").strip()
    name = taxon_element.findtext("ScientificName", default="").strip()

    if not rank or rank.lower() in EXCLUDED_RANKS:
        return None

    return {
        "rank": rank,
        "taxid": taxid,
        "name": name,
    }


def fetch_ncbi_lineage(taxid):
    """Fetch the ranked lineage for a taxon from NCBI Taxonomy."""
    query = urllib.parse.urlencode(
        {
            "db": "taxonomy",
            "id": taxid,
            "retmode": "xml",
        }
    )
    response = open_url(f"{NCBI_EFETCH_URL}?{query}", "application/xml")

    try:
        root = ET.fromstring(response)
    except ET.ParseError as error:
        raise RuntimeError("NCBI returned invalid taxonomy XML.") from error

    taxon = root.find("./Taxon")

    if taxon is None:
        raise RuntimeError(f"NCBI Taxonomy returned no record for taxid {taxid}.")

    lineage = []

    for ancestor in taxon.findall("./LineageEx/Taxon"):
        entry = parse_ncbi_taxon(ancestor)

        if entry:
            lineage.append(entry)

    current_taxon = parse_ncbi_taxon(taxon)

    if current_taxon:
        lineage.append(current_taxon)

    if not lineage:
        raise RuntimeError(f"No ranked lineage was found for taxid {taxid}.")

    return lineage


def write_lineage(lineage, output):
    """Write the taxonomy lineage as a tab-separated file."""
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "taxid", "name"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(lineage)


def fetch_uniprot_mnemonic(taxid):
    """Fetch the UniProt mnemonic for a taxon."""
    url = UNIPROT_TAXONOMY_URL.format(taxid=taxid)

    try:
        response = open_url(url, "application/json")
        data = json.loads(response)
    except (RuntimeError, json.JSONDecodeError) as error:
        print(
            f"Warning: failed to fetch UniProt mnemonic for taxid {taxid}: {error}",
            file=sys.stderr,
        )
        return "NA"

    return data.get("mnemonic") or "NA"


def main():
    args = parse_args()
    taxid = args.taxid.strip()

    if not taxid.isdigit():
        raise SystemExit("TaxID must contain only digits.")

    lineage_output = Path(f"{taxid}_lineage.tsv")
    mnemonic_output = Path(f"{taxid}.mnemonic")

    print(f"Fetching NCBI taxonomy lineage for taxid {taxid}")
    lineage = fetch_ncbi_lineage(taxid)
    write_lineage(lineage, lineage_output)
    print(f"Fetched {len(lineage)} ranked taxa")

    mnemonic = fetch_uniprot_mnemonic(taxid)
    mnemonic_output.write_text(f"{mnemonic}\n", encoding="utf-8")
    print(f"UniProt mnemonic: {mnemonic}")


if __name__ == "__main__":
    main()