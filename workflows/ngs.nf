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

include { PARSE_DESIGN_SWF                } from "../subworkflows/parse_design_SWF.nf"
include { PREPARE_REFS                    } from "../subworkflows/prepare_refs.nf"
include { cutadapt                        } from "../modules/cutadapt.nf"
include { fastqc as fastqc_raw            } from "../modules/fastqc.nf"
include { fastqc as fastqc_trim           } from "../modules/fastqc.nf"
include { multiqc as multiqc_raw          } from "../modules/multiqc.nf"
include { multiqc as multiqc_trim         } from "../modules/multiqc.nf"

workflow NGS {
    /*
    ---------------------------------------------------------------------
        Read design file, parse sample names and identifiers, and stage reads files
    ---------------------------------------------------------------------
    */

    // set channel for input design file
    ch_input = file(params.input)

    // Subworkflow: Parse design file
    PARSE_DESIGN_SWF(ch_input)
    ch_reads_raw = PARSE_DESIGN_SWF.out.samples


    /*
    ---------------------------------------------------------------------
        Prepare reference files
    ---------------------------------------------------------------------
    */
    PREPARE_REFS(
        params.ref,
        params.annot
    )
    ch_ref   = PREPARE_REFS.out.fasta
    ch_annot = PREPARE_REFS.out.annotations


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
