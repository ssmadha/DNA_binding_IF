#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DOWNLOAD_GENE {

    time { 15.m * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "results/individual", mode: 'copy'

    conda "environment.yaml"

    input:
    val gene_name

    output:
    tuple val(gene_name), path("${gene_name}.txt")

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

process LOG_FAILURES {

    publishDir "results", mode: 'copy'

    input:
    val failed_list

    output:
    path "failed_genes.txt"

    script:
    """
    printf "%s\n" ${failed_list.join(' ')} > failed_genes.txt
    """
}

workflow {

    genes = Channel
        .fromPath("expressed_coding_isoforms_with_relative_tpm_threshold_1_geneIDs.txt")
        .splitText()
        .map { it.trim() }
        .filter { it }

    gene_outputs = DOWNLOAD_GENE(genes)

    // Successful gene names
    successful_genes = gene_outputs.map { it[0] }

    // Collect lists
    all_genes_list = genes.collect()
    successful_list = successful_genes.collect()

    // Compute failures
    failed_genes = all_genes_list
        .combine(successful_list)
        .map { all, success -> all - success }

    COMBINE_GENES(gene_outputs.map{ it[1] }.collect())
    LOG_FAILURES(failed_genes)
}
