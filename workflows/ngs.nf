/*
=====================================================================
    NGS WORKFLOW
=====================================================================
*/

/*
---------------------------------------------------------------------
Import modules
---------------------------------------------------------------------
*/

include { CHECK_QUALITY                   } from "../workflows/check_quality.nf"
include { PREPARE_INPUTS                  } from "../workflows/prepare_inputs.nf"
include { PROCESS_READS                   } from "../workflows/process_reads.nf"

workflow NGS {
    PREPARE_INPUTS(
        params.samplesheet,
        params.genome,
        params.annotations
    )
    ch_reads_raw   = PREPARE_INPUTS.out.samples
    ch_genome      = PREPARE_INPUTS.out.genome
    ch_annotations = PREPARE_INPUTS.out.annotations

    PROCESS_READS(ch_reads_raw)
    ch_reads_pre_align = PROCESS_READS.out.reads_pre_align
    ch_trim_log        = PROCESS_READS.out.trim_log

    CHECK_QUALITY(
        ch_reads_raw,
        ch_reads_pre_align,
        ch_trim_log
    )

    /*
    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */

}
