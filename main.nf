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

import groovy.json.JsonSlurper

def jsonSlurper = new JsonSlurper()

LinkedHashMap defaults = jsonSlurper.parse(new File('/home/treevooor/SRAlign/parameter_specifications.json'))

println "Defaults"
println defaults

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
