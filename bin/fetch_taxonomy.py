#!/usr/bin/env python3
"""
Simple utility:
- fetch lineage from lbgi API for a taxid
- fetch mnemonic info from UniProt/EBI taxonomy API

Usage:
fetch_taxonomy.py <TAXID>
"""

import json
import sys
import urllib.request
import urllib.error
import csv

def fetch_lineage(taxid):
    """Return list of {rank,id,name} from lbgi lineage API or raise."""
    url = f"https://lbgi.fr/api/taxonomy/lineage/{taxid}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            data = json.load(r)
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Failed to fetch lineage: HTTP {e.code}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to fetch lineage: {e}") from e

    # Expect structure { "data": [ ... ] }
    lineage = []
    for item in data.get("data", []):
        rank = item.get("rank")
        # keep only entries with a rank and exclude 'clade' or 'no rank'
        if rank and rank.lower() != "clade" and rank.lower() != "no rank":
            lineage.append({"rank": rank, "id": item.get("id"), "name": item.get("name")})

    with open(f"{taxid}_lineage.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rank", "id", "name"], delimiter='\t')
        writer.writeheader()
        for item in lineage:
            # ensure values are strings
            writer.writerow({
                "rank": item.get("rank", ""),
                "id": item.get("id", ""),
                "name": item.get("name", "")
            })

    print(f"Fetched lineage with {len(lineage)} ranks.")
    return lineage

def fetch_mnemonic(taxid):
    """Fetch mnemonic from EBI/UniProt taxonomy API. Returns string or None.
    """
    url = f"https://www.ebi.ac.uk/proteins/api/taxonomy/id/{taxid}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            data = json.load(r)
    except Exception:
        print(f"Warning: failed to fetch mnemonic for taxid {taxid}", file=sys.stderr)
        return 'NA'

    mnemonic = data.get("mnemonic")
    print(f"Fetched mnemonic: {mnemonic}", file=sys.stderr)
    with open(f"{taxid}.mnemonic", "w", encoding="utf-8") as fh:
        fh.write(mnemonic if mnemonic else 'NA')

    return mnemonic


def main():
    taxid = sys.argv[1]
    print(f"Fetching info for taxid {taxid}")
    fetch_lineage(taxid)
    fetch_mnemonic(taxid)

if __name__ == "__main__":
    main()
