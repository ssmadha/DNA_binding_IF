#!/bin/bash

while read -r ensg_id ;
do
  python -m TF_ASIF.download_gene "$ensg_id"
done < TFs_Ensembl_v_1.01.txt