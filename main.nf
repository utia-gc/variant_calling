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
---------------------------------------------------------------------
    HELP MESSAGE
---------------------------------------------------------------------
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
---------------------------------------------------------------------
    SET GENOME PARAMETERS
---------------------------------------------------------------------
*/

def setGenomeParams (attribute) {
    return params.genomes[ params.genome ][ attribute ]
}


params.fasta   = setGenomeParams('fasta')
params.bowtie2 = setGenomeParams('bowtie2')
params.hisat2  = setGenomeParams('hisat2')

// set parameters for contamination check
def setContaminantParams (attribute) {
    return params.genomes[ params.contaminant ][ attribute ]
}


if (params.contaminant) {
    fastaContaminant   = setContamParams('fasta')
    bowtie2Contaminant = setContamParams('bowtie2')
    hisat2Contaminant  = setContamParams('hisat2')
}

/*
---------------------------------------------------------------------
    RUN MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { sralign } from './workflows/sralign.nf'

workflow {
    sralign()
}
