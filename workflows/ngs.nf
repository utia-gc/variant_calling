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
include { cutadapt                        } from "../modules/cutadapt.nf"
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


    /*
    ---------------------------------------------------------------------
        Reads QC
    ---------------------------------------------------------------------
    */

    if(!params.skipRawFastQC) {
        fastqc_raw(ch_reads_raw)
        multiqc_raw(fastqc_raw.out.zip.collect(), "raw")
    }

    /*
    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------
    */

    if(!params.skipTrimReads) {


        cutadapt(ch_reads_raw,
                          params.r1_adapter,
                          params.r2_adapter,
                          params.minimum_length)
        ch_reads_pre_align = cutadapt.out.reads

        if(!params.skipTrimFastQC) {
            fastqc_trim(ch_reads_pre_align)
            multiqc_trim(fastqc_trim.out.zip.collect(), "trimmed")
        }

    } else {
        ch_reads_pre_align = ch_reads_raw
    }


    /*
    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */

}
