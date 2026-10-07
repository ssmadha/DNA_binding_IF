#!/usr/bin/env python
"""
Merge per-gene download_gene.py outputs into one TSV.

Current per-gene or per-batch outputs (*.tsv) are already TSVs with a header (see
TF_ASIF.matching.COVERAGE_COLUMNS: one row per transcript and domain);
they are concatenated under a single header.

Older runs (ENSG*.txt) are download_gene.py's former printed output: a
header line with the gene ID, then one 4-line record per transcript as
printed by Transcript.align_to_reference:

    ENSG...                     gene ID
    ENST...                     transcript ID
    <int>                       number of superisoform domains
    [<float>, <float>, ...]     per-domain coverage of this transcript

These are parsed into columns gene_id, transcript_id, n_domains,
domain_coverage (domain_coverage kept in the same "[a, b, ...]" list form
as the input); transcripts with n_domains == 0 are dropped unless
--keep-empty is set.
"""
import argparse
import glob
import os
import re
import sys

GENE_RE = re.compile(r"^ENSG\d+$")
TRANSCRIPT_RE = re.compile(r"^ENST\d+$")
INT_RE = re.compile(r"^\d+$")
LIST_RE = re.compile(r"^\[.*\]$")

COLUMNS = ["gene_id", "transcript_id", "n_domains", "domain_coverage"]
# Header of current per-gene TSVs; same as TF_ASIF.matching.COVERAGE_COLUMNS
# (not imported, so this script keeps running without pandas/biopython).
COVERAGE_COLUMNS = ["gene_id", "transcript_id", "domain_id", "domain_type", "source_transcript_id",
                    "positions", "n_residues", "n_covered", "coverage"]


def parse_gene_file(path):
    """
    Parse one per-gene output file.

    Parameters
    ----------
    path: str
        Path to a download_gene.py output file.

    Returns
    -------
    tuple(list, list)
        (records, skipped): records is a list of
        (gene_id, transcript_id, n_domains, domain_coverage) tuples;
        skipped is a list of (line_number, line) that didn't fit the
        expected record layout.
    """
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle]

    records = []
    skipped = []
    # Line 1 is download_gene.py's own echo of the gene ID, not part of
    # a transcript record.
    i = 1 if lines and GENE_RE.match(lines[0]) else 0
    while i < len(lines):
        block = lines[i:i + 4]
        if (len(block) == 4 and GENE_RE.match(block[0]) and TRANSCRIPT_RE.match(block[1])
                and INT_RE.match(block[2]) and LIST_RE.match(block[3])):
            records.append((block[0], block[1], int(block[2]), block[3]))
            i += 4
        else:
            if lines[i].strip():
                skipped.append((i + 1, lines[i]))
            i += 1
    return records, skipped


def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input-dir", required=True,
                        help="Directory of per-gene/batch output files (*.tsv, or ENSG*.txt from older runs)")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to write the merged TSV to")
    parser.add_argument("--keep-empty", action="store_true",
                        help="Older (.txt) runs only: keep transcripts with n_domains == 0 (dropped by default)")
    return parser.parse_args()


def merge_tsv_files(gene_files, output):
    """
    Concatenate per-gene TSV outputs under one header.

    Parameters
    ----------
    gene_files: list[str]
    output: str
    """
    header = None
    n_rows = 0
    with open(output, "w") as out:
        for path in gene_files:
            with open(path) as handle:
                file_header = handle.readline()
                if header is None:
                    header = file_header
                    out.write(header)
                elif file_header != header:
                    sys.exit("Header of %s does not match the other files" % path)
                for line in handle:
                    if line.strip():
                        out.write(line)
                        n_rows += 1
    print("Merged %d gene files: wrote %d rows -> %s" % (len(gene_files), n_rows, output), file=sys.stderr)


if __name__ == "__main__":
    args = get_args()
    # Any .tsv, not just ENSG*.tsv: batch files are named after their
    # first gene list entry, which needn't be an Ensembl ID (e.g.
    # DUX1_HUMAN_batch.tsv); mismatched headers are still rejected.
    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*.tsv")))
    gene_files = sorted(glob.glob(os.path.join(args.input_dir, "ENSG*.txt")))
    if tsv_files and gene_files:
        sys.exit("Found both .tsv and older ENSG*.txt outputs in " + args.input_dir + "; merge them separately")
    if tsv_files:
        merge_tsv_files(tsv_files, args.output)
        sys.exit(0)
    if not gene_files:
        print("WARNING: no .tsv or ENSG*.txt files found in " + args.input_dir, file=sys.stderr)
        with open(args.output, "w") as out:
            out.write("\t".join(COVERAGE_COLUMNS) + "\n")
        sys.exit(0)

    n_written = 0
    n_empty = 0
    n_skipped = 0
    with open(args.output, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for path in gene_files:
            records, skipped = parse_gene_file(path)
            for line_number, line in skipped:
                print("WARNING: %s:%d: unexpected line skipped: %s" % (path, line_number, line),
                      file=sys.stderr)
            n_skipped += len(skipped)
            for gene_id, transcript_id, n_domains, domain_coverage in records:
                if n_domains == 0 and not args.keep_empty:
                    n_empty += 1
                    continue
                out.write("%s\t%s\t%d\t%s\n" % (gene_id, transcript_id, n_domains, domain_coverage))
                n_written += 1

    print("Merged %d gene files: wrote %d transcripts, dropped %d with no domains, "
          "skipped %d unexpected lines -> %s"
          % (len(gene_files), n_written, n_empty, n_skipped, args.output), file=sys.stderr)