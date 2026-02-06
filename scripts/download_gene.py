#!/home/shariq/anaconda3/envs/DNA_Binding_IF/bin/python

import sys

sys.path.append("/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/classes")
import gene

if __name__ == "__main__":
    argv = sys.argv
    gene = gene.Gene(argv[1])