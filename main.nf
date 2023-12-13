/*
---------------------------------------------------------------------
    utia-gc/ngs
---------------------------------------------------------------------
https://github.com/utia-gc/ngs
*/

nextflow.enable.dsl=2

/*
---------------------------------------------------------------------
    RUN MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { CHECK_QUALITY  } from "./workflows/check_quality.nf"
include { MAP_READS      } from "./workflows/map_reads.nf"
include { PREPARE_INPUTS } from "./workflows/prepare_inputs.nf"
include { PROCESS_READS  } from "./workflows/process_reads.nf"

/*
---------------------------------------------------------------------
    CHECK FOR REQUIRED PARAMETERS
---------------------------------------------------------------------
*/
if (params.samplesheet) {
    log.info "Using samplesheet '${params.samplesheet}'"
} else {
    log.error "Parameter 'samplesheet' is required but was not provided."
    exit 64
}
if (params.genome) {
    log.info "Using reference genome '${params.genome}'"
} else {
    log.error "Parameter 'genome' is required but was not provided."
    exit 64
}
if (params.annotations) {
    log.info "Using reference annotations '${params.annotations}'"
} else {
    log.error "Parameter 'annotations' is required but was not provided."
    exit 64
}

workflow {
    PREPARE_INPUTS(
        params.samplesheet,
        params.genome,
        params.annotations
    )
    ch_reads_raw    = PREPARE_INPUTS.out.samples
    ch_reads_raw.dump(tag: "ch_reads_raw")
    ch_genome       = PREPARE_INPUTS.out.genome
    ch_genome_index = PREPARE_INPUTS.out.genome_index
    ch_annotations  = PREPARE_INPUTS.out.annotations

    PROCESS_READS(ch_reads_raw)
    ch_reads_pre_align = PROCESS_READS.out.reads_pre_align
    ch_trim_log        = PROCESS_READS.out.trim_log

    MAP_READS(
        ch_reads_pre_align,
        ch_genome,
        ch_annotations,
        params.tools.map
    )
    ch_alignmentsIndividual = MAP_READS.out.alignmentsIndividual
    ch_alignmentsMerged     = MAP_READS.out.alignmentsMerged

    CHECK_QUALITY(
        ch_reads_raw,
        ch_reads_pre_align,
        ch_trim_log,
        ch_genome_index,
        ch_alignmentsIndividual,
        ch_alignmentsMerged
    )
}
