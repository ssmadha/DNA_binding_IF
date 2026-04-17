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
                            choices=["ppi", "dbi"],
                            default=["ppi"],
                            help="Which domain types to use. Options are %(choices)s. (Default: %(default)s)")

    return getoptions.parse_args()


if __name__ == "__main__":
    args = get_args()
    print(args.ensgid)
    if args.ensgid.startswith("ENSG"):
        test_gene = gene.Gene(args.ensgid, refmode=args.refmode, domain_filter=args.domains)