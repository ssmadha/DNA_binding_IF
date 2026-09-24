#!/usr/bin/env python
import argparse
import sys

import TF_ASIF.gene as gene

def get_args():
    getoptions = argparse.ArgumentParser()
    getoptions.add_argument("-e", "--ensgid",
                            required=True,
                            help="Ensembl Gene ID")
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
    getoptions.add_argument("-o", "--keepoverlappingdomains",
                            action="store_true",
                            help="Keep every classified domain from every transcript as-is, without collapsing \
                                  mutually overlapping domains (across all transcripts) down to the largest \
                                  domain in each overlapping group. Default: overlap collapsing is on.")

    return getoptions.parse_args()


if __name__ == "__main__":
    args = get_args()
    print(args.ensgid)
    if args.ensgid.startswith("ENSG"):
        test_gene = gene.Gene(args.ensgid, ppi_binding_site_file=args.ppibindingsitefile,
                              cds_fasta_file=args.cdsfastafile, uniprot_mapping_file=args.uniprotmappingfile,
                              gtf_file=args.gtffile, interpro_domains_file=args.interprodomainsfile,
                              refmode=args.refmode, domain_filter=args.domains,
                              merge_overlapping_domains=not args.keepoverlappingdomains)