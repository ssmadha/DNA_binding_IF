#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DOWNLOAD_GENE {

    time '15m'
    errorStrategy 'ignore'

    publishDir "results", mode: 'copy'

    conda "/home/shariq/anaconda3/envs/DNA_Binding_IF"

    input:
    val gene_name

    output:
    path "${gene_name}.txt"

    script:
    """
    download_gene.py ${gene_name} > ${gene_name}.txt
    """
}

process COMBINE_GENES {

    publishDir "results", mode: 'copy'

    input:
    path gene_files

    output:
    path "all_genes_combined.txt"

    script:
    """
    cat ${gene_files} > all_genes_combined.txt
    """
}

workflow {

    genes = Channel
        .fromPath("expressed_coding_isoforms_with_relative_tpm_threshold_1_geneIDs.txt")
        .splitText()
        .map { it.trim() }
        .filter { it }


    gene_outputs = DOWNLOAD_GENE(genes)

    COMBINE_GENES(gene_outputs.collect())
}
