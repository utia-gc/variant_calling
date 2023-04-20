/*
=====================================================================
    utia_rnaseq WORKFLOW
=====================================================================
*/

/*
---------------------------------------------------------------------
Import modules
---------------------------------------------------------------------
*/

include { PARSE_DESIGN_SWF                } from "../subworkflows/parse_design_SWF.nf"
include { PREPARE_REFS                    } from "../subworkflows/prepare_refs.nf"
include { CUTADAPT_ADAPTERS               } from "../modules/cutadapt.nf"
include { FASTQC as FQRAW                 } from "../modules/fastqc.nf"
include { FASTQC as FQTRIM                } from "../modules/fastqc.nf"
include { MULTIQC as MQRAW                } from "../modules/multiqc.nf"
include { MULTIQC as MQTRIM               } from "../modules/multiqc.nf"
include { STAR_INDEX ; STAR_MAP           } from "../modules/star.nf"
include { SAMTOOLS_SORT                   } from "../modules/samtools.nf"

workflow UTIA_RNASEQ {
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
        FQRAW(ch_reads_raw, "raw")
        MQRAW(FQRAW.out.fastq_ch.collect(), "raw")
    }

    /*
    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipTrimReads) {


        CUTADAPT_ADAPTERS(ch_reads_raw,
                          params.r1_adapter,
                          params.r2_adapter,
                          params.minimum_length)
        ch_reads_pre_align = CUTADAPT_ADAPTERS.out.reads

        if (!params.skipTrimFastQC) {
            FQTRIM(ch_reads_pre_align, "trimmed")
            MQTRIM(FQTRIM.out.fastq_ch.collect(), "trimmed")
        }

    } else {
        ch_reads_pre_align = ch_reads_raw
    }


    /*
    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */

    if (params.alignmentTool == "star") {
        STAR_INDEX(ch_ref, ch_annot)
        STAR_MAP(ch_reads_pre_align, STAR_INDEX.out.star_idx)
    }

    SAMTOOLS_SORT(STAR_MAP.out.star_bam)
    SAMTOOLS_SORT.out.sort_star_bam.view()


    // TODO: Add htseq-count to workflow

}
