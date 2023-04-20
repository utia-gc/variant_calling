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

include { PARSE_DESIGN_SWF                } from "${projectDir}/subworkflows/parse_design_SWF.nf"
include { CUTADAPT_ADAPTERS               } from "${projectDir}/modules/cutadapt.nf"
include { FASTQC as FQRAW                 } from "${projectDir}/modules/fastqc.nf"
include { FASTQC as FQTRIM                } from "${projectDir}/modules/fastqc.nf"
include { MULTIQC as MQRAW                } from "${projectDir}/modules/multiqc.nf"
include { MULTIQC as MQTRIM               } from "${projectDir}/modules/multiqc.nf"
include { STAR_INDEX ; STAR_MAP           } from "${projectDir}/modules/star.nf"
include { SAMTOOLS_SORT                   } from "${projectDir}/modules/samtools.nf"

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
    ch_reads_raw = PARSE_DESIGN_SWF.out.reads


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
        STAR_INDEX(params.ref, params.annot)
        STAR_MAP(ch_reads_pre_align, STAR_INDEX.out.star_idx)
    }

    SAMTOOLS_SORT(STAR_MAP.out.star_bam)
    SAMTOOLS_SORT.out.sort_star_bam.view()


    // TODO: Add htseq-count to workflow

}
