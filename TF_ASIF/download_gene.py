#!/home/shariq/anaconda3/envs/DNA_Binding_IF/bin/python

import sys

import TF_ASIF.gene as gene

if __name__ == "__main__":
    argv = sys.argv
    print(argv[1])
    if argv[1].startswith("ENSG"):
        test_gene = gene.Gene(argv[1])
        print(test_gene.superisoform_seq)