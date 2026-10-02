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

With --database pfam, pulls the raw Pfam hits instead (same schema,
source "Pfam") - used by Domain.determine_protein_interaction as
evidence for PPI domains, not as domains themselves.
"""
import argparse
import gzip
import csv

import requests

# Pinned to the Ensembl release 109 archive (Feb 2023), matching the
# release of the GTF/CDS/UniProt-xref files the pipeline also uses. The
# live www.ensembl.org mart tracks the current release instead.
BIOMART_URL = "https://feb2023.archive.ensembl.org/biomart/martservice"
DATASET = "hsapiens_gene_ensembl"
# BioMart attribute prefix, default output file and source label for
# each supported database.
DATABASES = {
    "interpro": ("interpro", "../../reference_data/Homo_sapiens.GRCh38.interpro_domains.tsv.gz", "SuperFamily"),
    "pfam": ("pfam", "../../reference_data/Homo_sapiens.GRCh38.pfam_domains.tsv.gz", "Pfam"),
}
OUTPUT_FILE = DATABASES["interpro"][1]

QUERY_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6" completionStamp="1">
  <Dataset name="{dataset}" interface="default">
    <Attribute name="ensembl_peptide_id"/>
    <Attribute name="{attribute}"/>
    <Attribute name="{attribute}_start"/>
    <Attribute name="{attribute}_end"/>
  </Dataset>
</Query>
"""

RAW_TSV_FILE = "interpro_domains_raw.tsv"


def download_raw_tsv(raw_tsv_file=RAW_TSV_FILE, biomart_url=BIOMART_URL, database="interpro"):
    query_xml = QUERY_TEMPLATE.format(dataset=DATASET, attribute=DATABASES[database][0])
    print("Querying BioMart at " + biomart_url + " for the whole human proteome's " + database + " domains "
          "(protein_id, domain_id, positions, source); this can take several minutes...")
    with requests.get(biomart_url, params={"query": query_xml}, stream=True, timeout=900) as response:
        response.raise_for_status()
        with open(raw_tsv_file, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                handle.write(chunk)
    # completionStamp="1" makes BioMart append a final "[success]" line;
    # without it, a server-side abort mid-export looks like a short but
    # otherwise valid HTTP 200 response.
    with open(raw_tsv_file, "rb") as handle:
        handle.seek(0, 2)
        handle.seek(max(0, handle.tell() - 64))
        tail = handle.read().decode(errors="replace").strip().splitlines()
    if not tail or tail[-1] != "[success]":
        raise RuntimeError("BioMart response is incomplete (no [success] stamp): " + raw_tsv_file)
    print("Download complete: " + raw_tsv_file)


def clean_and_compress(raw_tsv_file=RAW_TSV_FILE, output_file=OUTPUT_FILE, source="SuperFamily"):
    # Same column shape as reference_data/ppi_binding_sites.tsv
    # (protein_id, domain_id, positions, source) so both files can be
    # read by one shared parser. Each row here is always a single
    # contiguous span, so "positions" is just "start-end" - never a
    # multi-segment list like the PPI file's discontinuous contacts.
    seen = set()
    kept = 0
    total = 0
    with open(raw_tsv_file, newline="") as fin, gzip.open(output_file, "wt", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["protein_id", "domain_id", "positions", "source"])
        for row in csv.reader(fin, delimiter="\t"):
            if row == ["[success]"]:
                continue
            total += 1
            if len(row) != 4 or not all(row):
                continue
            key = tuple(row)
            if key in seen:
                continue
            seen.add(key)
            protein_id, domain_id, start, end = row
            writer.writerow([protein_id, domain_id, f"{start}-{end}", source])
            kept += 1
    print(f"Wrote {kept} rows (from {total} raw rows) to {output_file}")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", choices=sorted(DATABASES), default="interpro",
                        help="Which domain hits to download. Default: %(default)s")
    parser.add_argument("-o", "--output", default=None,
                        help="Path to write the gzipped domains TSV to. \
                              Default: the database's file in reference_data/")
    parser.add_argument("-r", "--raw-tsv-file", default=RAW_TSV_FILE,
                        help="Path to write the raw (pre-cleaning) BioMart TSV response to. \
                              Default: %(default)s")
    parser.add_argument("-u", "--biomart-url", default=BIOMART_URL,
                        help="BioMart martservice endpoint to query. Default: %(default)s")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    _, default_output, source = DATABASES[args.database]
    download_raw_tsv(raw_tsv_file=args.raw_tsv_file, biomart_url=args.biomart_url, database=args.database)
    clean_and_compress(raw_tsv_file=args.raw_tsv_file, output_file=args.output or default_output, source=source)