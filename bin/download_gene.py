#!/usr/bin/env python
"""
Score one gene, or a batch of genes, and write one coverage TSV (see
TF_ASIF.matching.COVERAGE_COLUMNS). Reference files are loaded once and
reused across the batch. A gene that raises an error or exceeds
--genetimeout is logged to stderr and skipped, so it doesn't take the
rest of the batch down; --donefile lists the genes that finished.
"""
import argparse
import csv
import signal
import sys
import traceback

import TF_ASIF.gene as gene
from TF_ASIF import matching
from TF_ASIF.matching import COVERAGE_COLUMNS
from TF_ASIF.transcript import Transcript

def get_args():
    getoptions = argparse.ArgumentParser()
    getoptions.add_argument("-e", "--ensgid",
                            nargs="+",
                            default=[],
                            help="Ensembl Gene ID(s)")
    getoptions.add_argument("-f", "--genesfile",
                            help="File of Ensembl Gene IDs, one per line (combined with any --ensgid)")
    getoptions.add_argument("--output",
                            help="Path to write the coverage TSV to. Default: stdout")
    getoptions.add_argument("--donefile",
                            help="Path to write the IDs of genes that finished (with or without rows) to")
    getoptions.add_argument("-t", "--genetimeout",
                            type=int,
                            default=0,
                            help="Skip a gene after this many seconds (0 = no limit). (Default: %(default)s)")
    getoptions.add_argument("-m", "--refmode",
                            default="superisoform",
                            choices=["superisoform", "reference"],
                            help="Whether to use a superisoform or a reference isoform for alignment. \
                                  Reference isoform not currently supported.")
    getoptions.add_argument("-d", "--domains",
                            nargs="*",
                            choices=["ppi_domain", "ppi_bs", "dbi"],
                            default=["ppi_domain", "dbi"],
                            help="Which domain types to use. Options are %(choices)s. (Default: %(default)s)")
    getoptions.add_argument("-b", "--ppibindingsitefile",
                            help="File with PPI binding sites")
    getoptions.add_argument("-c", "--cdsfastafile",
                            help="Ensembl CDS FASTA file, used for local protein sequence lookup")
    getoptions.add_argument("-u", "--uniprotmappingfile",
                            help="Ensembl protein-to-UniProt xref TSV file, used for UniProt ID resolution")
    getoptions.add_argument("-g", "--gtffile",
                            help="Ensembl GTF annotation file, used to list this gene's isoforms")
    getoptions.add_argument("-p", "--interprodomainsfile",
                            help="Local InterPro protein domain TSV file, used for SuperFamily/InterPro domain lookup")
    getoptions.add_argument("-P", "--pfamdomainsfile",
                            help="Local Pfam hit TSV file; hits of 3did families are the evidence for PPI domains. \
                                  Required with -d ppi_domain.")
    getoptions.add_argument("-T", "--interproentrytypesfile",
                            help="InterPro entry.list file (entry types); only Domain, Homologous_superfamily and \
                                  Repeat entries can be PPI domains. Required with -d ppi_domain.")
    getoptions.add_argument("-o", "--keepoverlappingdomains",
                            action="store_true",
                            help="Keep every classified domain from every transcript as-is, without collapsing \
                                  mutually overlapping domains (across all transcripts) down to the largest \
                                  domain in each overlapping group. Default: overlap collapsing is on.")

    getoptions.add_argument("-M", "--matchmode",
                            default="segment",
                            choices=["segment", "alignment"],
                            help="How transcripts are matched to domains: exact codon matching against a \
                                  segment superisoform, or the original protein alignment. (Default: %(default)s)")
    getoptions.add_argument("-i", "--identicalonly",
                            action="store_true",
                            help="Only count a reference residue as covered when it aligns to an identical \
                                  residue (mismatches count as uncovered, like gaps). Alignment mode only. \
                                  Default: any aligned residue counts as covered.")

    args = getoptions.parse_args()
    if not args.ensgid and not args.genesfile:
        getoptions.error("give at least one --ensgid or a --genesfile")
    if "ppi_domain" in args.domains and not (args.pfamdomainsfile and args.interproentrytypesfile):
        getoptions.error("-d ppi_domain needs --pfamdomainsfile and --interproentrytypesfile")
    return args


class GeneTimeout(Exception):
    pass


def _raise_timeout(signum, frame):
    raise GeneTimeout()


def preload_reference_files(args):
    """
    Load and cache the reference files this run's match mode reads, so
    that loading them isn't charged to the first gene's --genetimeout.
    """
    gene.Gene._get_gtf_index(args.gtffile)
    if "ppi_domain" in args.domains or "dbi" in args.domains:
        Transcript._get_interpro_domains(args.interprodomainsfile)
    if "ppi_domain" in args.domains:
        Transcript._get_ppi_regions(args.pfamdomainsfile)
        Transcript._get_interpro_entry_types(args.interproentrytypesfile)
    if "ppi_bs" in args.domains:
        Transcript._get_ppi_binding_sites(args.ppibindingsitefile)
    if args.matchmode == "segment":
        if "ppi_bs" in args.domains:
            matching.load_uniprot_xrefs(args.uniprotmappingfile)
    else:
        Transcript._get_uniprot_mapping(args.uniprotmappingfile)
        Transcript._get_cds_sequences(args.cdsfastafile)


if __name__ == "__main__":
    args = get_args()
    gene_ids = list(args.ensgid)
    if args.genesfile:
        with open(args.genesfile) as handle:
            gene_ids += [line.strip() for line in handle if line.strip()]
    signal.signal(signal.SIGALRM, _raise_timeout)
    preload_reference_files(args)

    # Coverage table goes to --output or stdout (TSV, header always
    # written); warnings and messages go to stderr.
    out = open(args.output, "w") if args.output else sys.stdout
    writer = csv.DictWriter(out, fieldnames=COVERAGE_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    done = []
    for gene_id in gene_ids:
        if not gene_id.startswith("ENSG"):
            print("Skipping %s: not an Ensembl Gene ID" % gene_id, file=sys.stderr)
            continue
        signal.alarm(args.genetimeout)
        try:
            test_gene = gene.Gene(gene_id, ppi_binding_site_file=args.ppibindingsitefile,
                                  cds_fasta_file=args.cdsfastafile, uniprot_mapping_file=args.uniprotmappingfile,
                                  gtf_file=args.gtffile, interpro_domains_file=args.interprodomainsfile,
                                  pfam_domains_file=args.pfamdomainsfile,
                                  interpro_entry_types_file=args.interproentrytypesfile,
                                  refmode=args.refmode, domain_filter=args.domains,
                                  merge_overlapping_domains=not args.keepoverlappingdomains,
                                  identical_only=args.identicalonly, matchmode=args.matchmode)
            signal.alarm(0)
        except GeneTimeout:
            print("FAILED %s: exceeded %d s" % (gene_id, args.genetimeout), file=sys.stderr)
            continue
        except Exception:
            signal.alarm(0)
            print("FAILED %s:\n%s" % (gene_id, traceback.format_exc()), file=sys.stderr)
            continue
        if not hasattr(test_gene, "strand"):
            print("FAILED %s: not found in %s" % (gene_id, args.gtffile), file=sys.stderr)
            continue
        writer.writerows(test_gene.coverage_rows)
        out.flush()
        done.append(gene_id)
    if args.output:
        out.close()
    if args.donefile:
        with open(args.donefile, "w") as handle:
            handle.writelines(gene_id + "\n" for gene_id in done)
    print("Finished %d of %d genes" % (len(done), len(gene_ids)), file=sys.stderr)