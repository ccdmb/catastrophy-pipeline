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


def handle_dbcan_params(params) {
    if ( (params.dbcan || params.dbcan_url) && !params.dbcan_version ) {
        exit 1, "Please provide the dbcan version that you are providing " +
                "using '--dbcan_version'."
    }

    if ( params.dbcan_version && !(params.dbcan_version in DBCAN_VERSIONS) ) {
        exit 1, "The dbcan version you provided is not supported. " +
                "Valid options are ${DBCAN_VERSIONS}."
    }

    // Set default version
    dbcan_version = params.dbcan_version ?: DBCAN_VERSIONS[-1]

    if ( params.dbcan ) {
        dbcan_file = file(params.dbcan, checkIfExists: true)
    } else {
        url = params.dbcan_url ?: DBCAN_URLS[dbcan_version]
        dbcan_file = download_dbcan(url)
    }

    return [dbcan_version, dbcan_file]
}


def handle_nomenclature_params(params) {
    if ( !(params.nomenclature in NOMENCLATURES) ) {
        exit 1, "The nomenclature you selected is not supported. " +
                "Valid options are ${NOMENCLATURES}."
    }

    nomenclature = params.nomenclature

    return nomenclature
}


def handle_proteomes_params(params) {

    if (params.proteomes) {
        // Can't use the .set {} syntax inside functions.
        ch_proteomes = Channel
            .fromPath(params.proteomes, checkIfExists: true, type: 'file')
            .map { f -> [f.simpleName, f] }
    } else {
        exit 1, "Please provide one or more fasta files to the --proteome parameter."
    }

    return ch_proteomes
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


process hmmpress {

    label "hmmer"
    label "small_task"

    input:
    path "database.txt"

    output:
    tuple path("database.txt"),
          path("database.txt.h3f"),
          path("database.txt.h3i"),
          path("database.txt.h3m"),
          path("database.txt.h3p"), emit: pressed_hmms

    script:
    """
    hmmpress database.txt
    """
}


process hmmscan {

    label "hmmer"
    label "small_task"

    tag "${name}"

    input:
    tuple path("database.txt"),
          path("database.txt.h3f"),
          path("database.txt.h3i"),
          path("database.txt.h3m"),
          path("database.txt.h3p")
    tuple val(name), path(proteome)

    output:
    tuple val(name), path("${name}.csv"), emit: domtab
    tuple val(name), path("${name}.txt"), emit: txt

    script:
    """
    hmmscan \
        --domtblout "${name}.csv" \
        database.txt \
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


workflow search_proteomes {

    get:
    database_file
    proteomes_ch

    main:
    pressed_hmms = hmmpress(database_file)
    (hmmscan_domtabs, hmmscan_txts) = hmmscan(pressed_hmms, proteomes_ch)

    emit:
    hmmscan_domtabs
    hmmscan_txts
}


workflow classify_proteomes {

    get:
    dbcan_version
    nomenclature
    domtabs

    main:
    transposed_domtabs = domtabs
        .toList()
        .map { it.transpose() }

    classified = catastrophy(dbcan_version, nomenclature, transposed_domtabs)

    emit:
    classified
}


workflow {

    main:
    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    proteomes_ch = handle_proteomes_params(params)
    (dbcan_version, dbcan_file) = handle_dbcan_params(params)
    nomenclature = handle_nomenclature_params(params)

    (hmmscan_domtabs, hmmscan_txts) = search_proteomes(dbcan_file, proteomes_ch)
    //classifications_file = classify_proteomes(dbcan_version, nomenclature, hmmscan_domtabs)

    publish:
    hmmscan_domtabs to: "${params.outdir}/matches"
    hmmscan_txts to: "${params.outdir}/matches"
    //classifications_file to: "${params.outdir}"
}
