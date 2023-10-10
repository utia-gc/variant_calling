include { Parse_Samplesheet } from "../subworkflows/parse_samplesheet.nf"
include { Prepare_Refs } from "../subworkflows/prepare_refs.nf"

/**
 * Workflow to handle and prepare input files.
 * This includes parsing samplesheet, decompressing archives, standardizing references, etc.
 * Essentially, if a file is passed to a process downstream in the pipeline, it should run through this workflow.
 */

workflow PREPARE_INPUTS {
    take:
        samplesheet
        genome
        annotations

    main:
        Parse_Samplesheet(samplesheet)

        Prepare_Refs(
            genome,
            annotations
        )
    
    emit:
        samples      = Parse_Samplesheet.out.samples
        genome       = Prepare_Refs.out.genome
        genome_index = Prepare_Refs.out.genome_index
        annotations  = Prepare_Refs.out.annotations
}
