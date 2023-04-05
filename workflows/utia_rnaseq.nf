/*
=====================================================================
    utia_rnaseq WORKFLOW
=====================================================================
*/

/*
This object takes care of many necessary steps upon construction:
    - Logs a header for the pipeline that prints pipeline name and logo
    - Prints a help message if help parameter is specified
    - Checks parameters
*/ 
//def srawf = new SRAlignWorkflow(log, params, workflow)


/*
---------------------------------------------------------------------
Import modules
---------------------------------------------------------------------
*/

include { ParseDesignSWF as ParseDesign } from "${projectDir}/subworkflows/ParseDesignSWF.nf"


workflow UTIA_RNASEQ {
    /*
    ---------------------------------------------------------------------
        Read design file, parse sample names and identifiers, and stage reads files
    ---------------------------------------------------------------------
    */

    // set channel for input design file
    ch_input = file(params.input)

    // Subworkflow: Parse design file
    ParseDesign(ch_input)
    ch_readsRaw = ParseDesign.out.reads
    ch_readsRaw.view()


    /*
    ---------------------------------------------------------------------
        Reads QC
    ---------------------------------------------------------------------

    if (!params.skipReadsQC) {
        // perform QC
    } else {
        // don't perform QC
    }


    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------

    if (!params.skipTrimReads) {
        // perform trimming and adapter removal
    } else {
        // don't perform trimming and adapter removal
    }


    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */
}
