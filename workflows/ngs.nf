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
include { CUTADAPT                        } from "../modules/cutadapt.nf"
include { FASTQC as FQRAW                 } from "../modules/fastqc.nf"
include { FASTQC as FQTRIM                } from "../modules/fastqc.nf"
include { MULTIQC as MQRAW                } from "../modules/multiqc.nf"
include { MULTIQC as MQTRIM               } from "../modules/multiqc.nf"

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

    if (!params.skipRawFastQC) {
        FQRAW(ch_reads_raw)
        MQRAW(FQRAW.out.zip.collect(), "raw")
    }

    /*
    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipTrimReads) {


        CUTADAPT(ch_reads_raw,
                          params.r1_adapter,
                          params.r2_adapter,
                          params.minimum_length)
        ch_reads_pre_align = CUTADAPT.out.reads

        if (!params.skipTrimFastQC) {
            FQTRIM(ch_reads_pre_align)
            MQTRIM(FQTRIM.out.zip.collect(), "trimmed")
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
