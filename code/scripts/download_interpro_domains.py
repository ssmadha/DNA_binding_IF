"""
Download and process the local InterPro protein domain reference file
(../../reference_data/Homo_sapiens.GRCh38.interpro_domains.tsv.gz) used
by Transcript.download_domains in place of the Ensembl REST
overlap/translation endpoint.

Pulls Ensembl's own consolidated InterPro domain table (the same
underlying data the REST API's "interpro" field is sourced from) via a
single BioMart attribute query. Note: BioMart's separate "superfamily"
attributes (superfamily/superfamily_start/superfamily_end) cannot be
combined with "interpro" attributes in one query - they come from
different underlying tables, and BioMart silently returns their full
cross product instead of a real per-hit join. Querying "interpro"
attributes alone avoids that trap.
"""
import argparse
import gzip
import csv

import requests

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"
DATASET = "hsapiens_gene_ensembl"
OUTPUT_FILE = "../../reference_data/Homo_sapiens.GRCh38.interpro_domains.tsv.gz"

QUERY_XML = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="{DATASET}" interface="default">
    <Attribute name="ensembl_peptide_id"/>
    <Attribute name="interpro"/>
    <Attribute name="interpro_start"/>
    <Attribute name="interpro_end"/>
  </Dataset>
</Query>
"""

RAW_TSV_FILE = "interpro_domains_raw.tsv"


def download_raw_tsv(raw_tsv_file=RAW_TSV_FILE):
    print("Querying BioMart for the whole human proteome's InterPro domains "
          "(protein_stable_id, interpro_id, start, end); this can take several minutes...")
    with requests.get(BIOMART_URL, params={"query": QUERY_XML}, stream=True, timeout=900) as response:
        response.raise_for_status()
        with open(raw_tsv_file, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                handle.write(chunk)
    print("Download complete: " + raw_tsv_file)


def clean_and_compress(raw_tsv_file=RAW_TSV_FILE, output_file=OUTPUT_FILE):
    seen = set()
    kept = 0
    total = 0
    with open(raw_tsv_file, newline="") as fin, gzip.open(output_file, "wt", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["protein_stable_id", "interpro_id", "start", "end"])
        for row in csv.reader(fin, delimiter="\t"):
            total += 1
            if len(row) != 4 or not all(row):
                continue
            key = tuple(row)
            if key in seen:
                continue
            seen.add(key)
            writer.writerow(row)
            kept += 1
    print(f"Wrote {kept} rows (from {total} raw rows) to {output_file}")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default=OUTPUT_FILE,
                        help="Path to write the gzipped InterPro domains TSV to. \
                              Default: %(default)s")
    parser.add_argument("-r", "--raw-tsv-file", default=RAW_TSV_FILE,
                        help="Path to write the raw (pre-cleaning) BioMart TSV response to. \
                              Default: %(default)s")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    download_raw_tsv(raw_tsv_file=args.raw_tsv_file)
    clean_and_compress(raw_tsv_file=args.raw_tsv_file, output_file=args.output)