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

include { PREPARE_INPUTS                  } from "../workflows/prepare_inputs.nf"
include { PROCESS_READS                   } from "../workflows/process_reads.nf"
include { fastqc as fastqc_raw            } from "../modules/fastqc.nf"
include { fastqc as fastqc_trim           } from "../modules/fastqc.nf"
include { multiqc as multiqc_raw          } from "../modules/multiqc.nf"
include { multiqc as multiqc_trim         } from "../modules/multiqc.nf"

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

    /*
    ---------------------------------------------------------------------
        Reads QC
    ---------------------------------------------------------------------
    */

    if(!params.skipRawFastQC) {
        fastqc_raw(ch_reads_raw)
        multiqc_raw(fastqc_raw.out.zip.collect(), "raw")
    }

    if(!params.skipTrimReads) {
        if(!params.skipTrimFastQC) {
            fastqc_trim(ch_reads_pre_align)
            multiqc_trim(
                Channel.empty()
                    .concat(fastqc_trim.out.zip)
                    .concat(ch_trim_log)
                    .collect(),
                "trimmed"
            )
        }
    }

    /*
    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */

}
