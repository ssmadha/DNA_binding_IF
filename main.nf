#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DOWNLOAD_GTF {

    storeDir "${projectDir}/reference_data"

    output:
    path "Homo_sapiens.GRCh38.109.gtf.gz"

    script:
    """
    curl -fsSL -o Homo_sapiens.GRCh38.109.gtf.gz https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
    """
}

process DOWNLOAD_CDS_FASTA {

    storeDir "${projectDir}/reference_data"

    output:
    path "Homo_sapiens.GRCh38.cds.all.fa.gz"

    script:
    """
    curl -fsSL -o Homo_sapiens.GRCh38.cds.all.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
    """
}

process DOWNLOAD_UNIPROT_MAPPING {

    storeDir "${projectDir}/reference_data"

    output:
    path "Homo_sapiens.GRCh38.109.uniprot.tsv.gz"

    script:
    """
    curl -fsSL -o Homo_sapiens.GRCh38.109.uniprot.tsv.gz https://ftp.ensembl.org/pub/release-109/tsv/homo_sapiens/Homo_sapiens.GRCh38.109.uniprot.tsv.gz
    """
}

process DOWNLOAD_INTERPRO_DOMAINS {

    // No static FTP file for this one - it's built from a live BioMart
    // query (code/scripts/download_interpro_domains.py), which takes
    // several minutes and can be flaky over that long a request.
    time { 30.m * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    storeDir "${projectDir}/reference_data"

    conda "environment.yaml"

    output:
    path "Homo_sapiens.GRCh38.interpro_domains.tsv.gz"

    script:
    """
    python3 ${projectDir}/code/scripts/download_interpro_domains.py --output Homo_sapiens.GRCh38.interpro_domains.tsv.gz
    """
}

process DOWNLOAD_GENE {

    time { 15.m * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "results/individual", mode: 'copy'

    conda "environment.yaml"

    input:
    val gene_name
    // Not referenced by name in the script below (it uses the
    // params.*_file paths directly), but declaring them as inputs
    // here makes this process wait on DOWNLOAD_GTF/DOWNLOAD_CDS_FASTA/
    // DOWNLOAD_UNIPROT_MAPPING/DOWNLOAD_INTERPRO_DOMAINS actually
    // finishing (or cache-hitting via storeDir) before it runs, instead
    // of racing them.
    path gtf_file
    path cds_fasta_file
    path uniprot_mapping_file
    path interpro_domains_file

    output:
    tuple val(gene_name), path("${gene_name}.txt")

    script:
    """
    download_gene.py -e ${gene_name} -b ${params.binding_site_file} -c ${params.cds_fasta_file} -u ${params.uniprot_mapping_file} -g ${params.gtf_file} -p ${params.interpro_domains_file}> ${gene_name}.txt
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
        .fromPath(params.genes_file)
        .splitText()
        .map { it.trim() }
        .filter { it }

    gtf_file = DOWNLOAD_GTF()
    cds_fasta_file = DOWNLOAD_CDS_FASTA()
    uniprot_mapping_file = DOWNLOAD_UNIPROT_MAPPING()
    interpro_domains_file = DOWNLOAD_INTERPRO_DOMAINS()

    gene_outputs = DOWNLOAD_GENE(genes, gtf_file, cds_fasta_file, uniprot_mapping_file, interpro_domains_file)

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