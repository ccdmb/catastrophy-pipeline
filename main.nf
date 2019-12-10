#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// Set some globals
DBCAN_VERSIONS = ["v6", "v7", "v8"]
DBCAN_URLS = [
    "v6": "http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V6.txt",
    "v7": "http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V7.txt",
    "v8": "http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt",
]

NOMENCLATURES = ["nomenclature1", "nomenclature2", "nomenclature3"]

def helpMessage() {
    log.info "# CATAStrophy-pipeline"
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/catastroflow --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --proteomes                   Path to input proteomes as fasta files (Glob patterns must be surrounded with quotes).
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}


process download_dbcan {
    label "download"
    label "small_task"

    input:
    val url

    output:
    path 'dbcan.txt', emit: downloaded_db

    script:
    """
    wget -O dbcan.txt "${url}"
    """
}

process press_hmms {
    label "hmmer"
    label "small_task"

    input:
    path "dbcan.txt"

    output:
    tuple path("dbcan.txt"),
          path("dbcan.txt.h3f"),
          path("dbcan.txt.h3i"),
          path("dbcan.txt.h3m"),
          path("dbcan.txt.h3p"), emit: pressed_hmms

    script:
    """
    hmmpress dbcan.txt
    """
}

process hmmscan {
    label "hmmer"
    label "small_task"

    tag "${name}"

    input:
    tuple path("dbcan.txt"),
          path("dbcan.txt.h3f"),
          path("dbcan.txt.h3i"),
          path("dbcan.txt.h3m"),
          path("dbcan.txt.h3p")
    tuple val(name), path(proteome)

    output:
    tuple val(name), path("${name}.csv"), emit: domtab
    tuple val(name), path("${name}.txt"), emit: txt

    script:
    """
    hmmscan \
        --domtblout "${name}.csv" \
        dbcan.txt \
        "${proteome}" \
    > "${name}.txt"
    """
}


process catastrophy {
    label "catastrophy"
    label "small_task"

    input:
    val version
    val nomenclature
    tuple val(names), path(domtabs)

    output:
    path "catastrophy.tsv"

    script:
    """
    catastrophy \
        --format "hmmer_domtab" \
        --model "${version}" \
        --nomenclature "${nomenclature}" \
        --outfile "catastrophy.tsv" \
        --label ${names} \
        -- \
        ${domtabs}
    """
}


workflow {

    main:
        // Show help message
        if (params.help) {
            helpMessage()
            exit 0
        }

        // Validate input and setup channels.
        if (params.proteomes) {
            Channel
                .fromPath(params.proteomes, checkIfExists: true, type: 'file')
                .map { f -> [f.simpleName, f] }
                .set { ch_proteomes }
        } else {
            exit 1, "Please provide one or more fasta files to the --proteome parameter."
        }

        if ( (params.dbcan || params.dbcan_url) && !params.dbcan_version ) {
            exit 1, "Please provide the dbcan version that you are providing " +
                    "using '--dbcan_version'."
        }

        if ( !(params.dbcan_version in DBCAN_VERSIONS) ) {
            exit 1, "The dbcan version you provided is not supported. " +
                    "Valid options are ${DBCAN_VERSIONS}."
        }

        if ( !(params.nomenclature in NOMENCLATURES) ) {
            exit 1, "The nomenclature you selected is not supported. " +
                    "Valid options are ${NOMENCLATURES}."
        }

        // Set default version
        dbcan_version = params.dbcan_version ?: DBCAN_VERSIONS[-1]
        nomenclature = params.nomenclature

        if ( params.dbcan ) {
            dbcan_file = file(params.dbcan, checkIfExists: true)
        } else {
            url = params.dbcan_url ?: DBCAN_URLS[dbcan_version]
            dbcan_file = download_dbcan(url)
        }

        pressed_hmms = press_hmms(dbcan_file)
        (hmmscan_domtabs, hmmscan_txts) = hmmscan(pressed_hmms, ch_proteomes)
        transposed_domtabs = hmmscan_domtabs
            .toList()
            .map { it.transpose() }
            .first()
        classified = catastrophy(dbcan_version, nomenclature, transposed_domtabs)

    publish:
        hmmscan_domtabs to: "${params.outdir}/matches"
        hmmscan_txts to: "${params.outdir}/matches"
        classified to: "${params.outdir}"
}
