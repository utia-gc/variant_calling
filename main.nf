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
    // I/O
    input         : [
        description : "Path to input csv file.",
        required    : true
    ],
    baseDirGenome : [
        default     : "${launchDir}/data/references",
        description : "Path to base directory where genomes are found."
    ],
    baseDirReport : [
        default     : "${launchDir}/results/reports",
        description : "Path to base directory where output reports will be stored."
    ],
    baseDirData   : [
        default     : "${launchDir}/results/data",
        description : "Path to base directory where output data will be stored."
    ],
    multiqcConfig : [
        default     : "${projectDir}/configs/full_multiqc_config.yaml",
        description : "Path to MultiQC config file."
    ],

    // reference
    genome        : [
        description : "Reference genome ID.",
        options     : ["TODO"],
        required    : true
    ],
    contaminant   : [
        default     : "EB1",
        description : "Contaminant genome ID.",
        options     : ["TODO"]
    ],

    // trim reads 
    trimTool      : [
        default     : "fastp",
        description : "Tool to use for read trimming.",
        options     : ["TODO"]
    ],
    adapterR1     : [
        default     : "",
        description : "Adapter sequence to trim from read 1."
    ],
    adapterR2     : [
        default     : "",
        description : "Adapter sequence to trim from read 2."
    ],
    skipTrimReads : [
        default     : false,
        description : "Skip the read trimming step.",
        skip        : true
    ],
    skipTrimReadsQC : [
        default     : false,
        description : "Skip QC report for trimmed reads.",
        skip        : true
    ],

    // alignment
    alignmentTool : [
        default     : "bowtie2",
        description : "Tool to use for alignment to reference. The same tool is used to align to both genome and contaminant."
    ],
    forceAlignRawReads : [
        default     : false,
        description : "Force alignment of raw reads even if read trimming is performed."
    ],
    skipAlignGenome : [
        default     : false,
        description : "Skip alignment of reads to reference.",
        skip        : true
    ],
    skipSamStatsQC : [
        default     : false,
        description : "Skip QC report for alignment.",
        skip        : true
    ],
    skipAlignContam : [
        default     : false,
        description : "Skip alignment of reads to contaminant.",
        skip        : true
    ],

    // library complexity
    skipPreseq    : [
        default     : false,
        description : "Skip Preseq.",
        skip        : true
    ],

    // miscellaneous params
    help          : [
        default     : false,
        description : "Display this help message and exit the program."
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
