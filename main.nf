#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Domain types passed to download_gene.py -d (space-separated; any of
// ppi_domain, ppi_bs, dbi). e.g. --domains ppi_bs for PPI binding sites only.
// Where results are published; set per run (e.g. --outdir results_ppi)
// to avoid overwriting earlier runs.
params.outdir = "results"
params.domains = "ppi_domain dbi"
// true passes --keepoverlappingdomains (skip collapsing overlapping domains).
params.keep_overlapping_domains = false
// true passes --identicalonly (count only identical aligned residues as covered).
params.identical_only = false
// "segment" (exact codon matching) or "alignment" (original protein alignment).
params.match_mode = "segment"
// Comma-separated transcript expression TSV(s) (gene_id, transcript_id,
// <tissue>_TPM, ...); one ASIF table is written per file. Unset skips ASIF.
params.expression_files = null
// ASIF impact factor = 1 - mean(sigmoid(alpha * (coverage - beta))).
params.asif_alpha = 63
params.asif_beta = 0.3

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

    conda "${projectDir}/environment.yaml"

    output:
    path "Homo_sapiens.GRCh38.interpro_domains.tsv.gz"

    script:
    """
    python3 ${projectDir}/code/scripts/download_interpro_domains.py --output Homo_sapiens.GRCh38.interpro_domains.tsv.gz
    """
}

process DOWNLOAD_GENE {

    time { 15.m * task.attempt }
    // Retry twice, then skip the gene rather than aborting the whole run;
    // skipped genes end up in <outdir>/failed_genes.txt.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    publishDir "${params.outdir}/individual", mode: 'copy'

    conda "${projectDir}/environment.yaml"

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
    tuple val(gene_name), path("${gene_name}.tsv")

    script:
    """
    download_gene.py ${params.keep_overlapping_domains ? '--keepoverlappingdomains' : ''} ${params.identical_only ? '--identicalonly' : ''} --matchmode ${params.match_mode} -d ${params.domains} -e ${gene_name} -b ${params.ppi_binding_site_file} -c ${params.cds_fasta_file} -u ${params.uniprot_mapping_file} -g ${params.gtf_file} -p ${params.interpro_domains_file} > ${gene_name}.tsv
    """
}

process MERGE_RESULTS {

    publishDir params.outdir, mode: 'copy'

    input:
    path gene_files, stageAs: "individual/*"

    output:
    path "all_results.tsv"

    // Reads the staged directory rather than taking the files as
    // arguments, so tens of thousands of genes don't hit ARG_MAX.
    script:
    """
    merge_gene_results.py --input-dir individual --output all_results.tsv
    """
}

process COMPUTE_ASIF {

    publishDir params.outdir, mode: 'copy'

    conda "${projectDir}/environment.yaml"

    input:
    path all_results
    path expression_file

    output:
    path "${expression_file.baseName}_ASIF.tsv"

    script:
    """
    compute_asif.py --results ${all_results} --expression ${expression_file} --alpha ${params.asif_alpha} --beta ${params.asif_beta} --output ${expression_file.baseName}_ASIF.tsv
    """
}

process LOG_FAILURES {

    publishDir params.outdir, mode: 'copy'

    input:
    val failed_list

    output:
    path "failed_genes.txt"

    script:
    if (failed_list)
        """
        printf "%s\n" ${failed_list.join(' ')} > failed_genes.txt
        """
    else
        """
        touch failed_genes.txt
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

    // Genes with no output after all retries. ifEmpty([]) keeps these
    // steps running even if every gene fails; wrapping each list in [ ]
    // stops combine() from flattening the two lists into one.
    all_genes_list = genes.collect().map { [it] }
    successful_list = gene_outputs.map { it[0] }.collect().ifEmpty([]).map { [it] }
    failed_genes = all_genes_list
        .combine(successful_list)
        .map { all, success -> all - success }

    all_results = MERGE_RESULTS(gene_outputs.map{ it[1] }.collect().ifEmpty([]))

    if (params.expression_files) {
        expression_files = Channel.fromPath(params.expression_files.tokenize(','), checkIfExists: true)
        COMPUTE_ASIF(all_results, expression_files)
    }
    LOG_FAILURES(failed_genes)
}
