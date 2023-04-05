/*
---------------------------------------------------------------------
    trev-f/utia_rnaseq
---------------------------------------------------------------------
https://github.com/trev-f/utia_rnaseq
*/

nextflow.enable.dsl=2

/*
---------------------------------------------------------------------
    RUN MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { UTIA_RNASEQ } from './workflows/utia_rnaseq.nf'

workflow {
    UTIA_RNASEQ()
}
