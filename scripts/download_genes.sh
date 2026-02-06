#!/bin/bash

cd ../

while read ensg_id ;
do
  ./scripts/download_gene.py $ensg_id
  break
done < TFs_Ensembl_v_1.01.txt