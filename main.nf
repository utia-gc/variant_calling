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
    HANDLE INPUTS
---------------------------------------------------------------------
*/
