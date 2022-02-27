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

/*
    Design and Inputs
*/

// check design file
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input design file not specified!'
}


// set input design name
inName = params.input.take(params.input.lastIndexOf('.')).split('/')[-1]


/*
    Genomes and References
*/

// check genome
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "Genome '${params.genome}' is not available.\n\tAvailable genomes include: ${params.genomes.keySet().join(", ")}"
}


// set genome variables
bt2Index = params.genome ? params.genomes[params.genome].bowtie2 ?: false : false


/*
---------------------------------------------------------------------
    MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { ParseDesignSWF  as ParseDesign  } from './subworkflows/ParseDesignSWF.nf'
include { RawReadsQCSWF   as RawReadsQC   } from './subworkflows/RawReadsQCSWF.nf'
include { TrimReadsSWF    as TrimReads    } from './subworkflows/TrimReadsSWF.nf'
include { TrimReadsQCSWF  as TrimReadsQC  } from './subworkflows/TrimReadsQCSWF.nf'
include { AlignBowtie2SWF as AlignBowtie2 } from './subworkflows/AlignBowtie2SWF.nf'


workflow {
    ParseDesign(ch_input)
    RawReadsQC(ParseDesign.out.rawReads, inName)
    TrimReads(ParseDesign.out.rawReads)
    TrimReadsQC(TrimReads.out.trimReads, inName) 

    AlignBowtie2(TrimReads.out.trimReads, inName, bt2Index)
}
