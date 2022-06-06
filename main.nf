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
    CREATE SRALIGNWORKFLOW OBJECT
---------------------------------------------------------------------
*/

/*
This object takes care of many necessary steps upon construction:
    - Logs a header for the pipeline that prints pipeline name and logo
    - Prints a help message if help parameter is specified
*/ 

LinkedHashMap defaults = [
    // required param
    genome        : [
        default     : "WBcel235",
        required    : true,
        description : "Identifier for reference genome. Must be an acceptable genome key.",
        options     : ["WBcel235", "EB1"]
    ],
    // optional param
    contaminant   : [
        default     : "EB1",
        description : "Identifier for contaminant genome. Must be an acceptable genome key.",
        options     : ["WBcel235", "EB1"]
    ],
    // skip param
    skipTrimReads : [
        default     : false,
        description : "Skip the read trimming step.",
        options     : [],
        skip        : true
    ]
]

def srawf = new SRAlignWorkflow(log, params, defaults)


/*
---------------------------------------------------------------------
    RUN MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { sralign } from './workflows/sralign.nf'

workflow {
    sralign()
}
