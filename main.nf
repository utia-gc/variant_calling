#!/usr/bin/env nextflow

/*
---------------------------------------------------------------------
trev-f/sralign
---------------------------------------------------------------------
sralign - A flexible pipeline for short read alignment to a reference.
https://github.com/trev-f/sralign
*/

nextflow.enable.dsl=2

/*
Help message
*/

// Write help message

def help_message() {
    log.info"""
    Usage:

        nextflow run sralign -profile docker --input YYYYMMDD_input.csv --genome WBCel235

    Required arguments:
        -profile
        --input
        --genome
    """.stripIndent()
}

// Show help message

if (params.help) {
    help_message()
    exit 0
}

/*
Handle input
*/

channel
    .fromPath(params.input)
    .splitCsv(header: false)
    .map { row -> ["${row[0]}_rep${row[1]}", [file(row[2]), file(row[3])]] }
    .set { ch_fastqQC }

/*
fastQC
*/

process fastQC {
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(name), file(fastq_file)

    output:
    publishDir "${PWD}/fastQC", mode: 'copy'
    path "*", emit: ch_fastQC

    script:
    """
    fastqc ${fastq_file}
    """
}

process multiQC {
    container 'ewels/multiqc:v1.11'

    input:
    path fastqc

    output:
    publishDir "${PWD}/fastQC", mode: 'copy'
    path "*"

    script:
    """
    multiqc ${fastqc}
    """
}

workflow {
    fastQC(ch_fastqQC)
    multiQC(fastQC.out.collect())
}
