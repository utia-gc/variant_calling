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
include { FASTQC                        } from "${projectDir}/modules/fastqc.nf"
include { MULTIQC                       } from "${projectDir}/modules/multiqc.nf"
include { CUTADAPT_ADAPTERS             } from "${projectDir}/modules/cutadapt.nf"
include { STAR_INDEX ; STAR_MAP         } from "${projectDir}/modules/star.nf"
include { SAMTOOLS_SORT                 } from "${projectDir}/modules/samtools.nf"

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


    /*
    ---------------------------------------------------------------------
        Reads QC
    ---------------------------------------------------------------------
    */

    if (!params.skipRawFastQC) {
        FASTQC(ch_readsRaw, "raw")
        MULTIQC(FASTQC.out.fastq_ch.collect(), "raw")
    }

    /*
    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipTrimReads) {
        CUTADAPT_ADAPTERS(ch_readsRaw,
                          params.r1_adapter,
                          params.r2_adapter,
                          params.minimum_length)
        ch_reads_pre_align = CUTADAPT_ADAPTERS.out.reads

        if (!params.skipTrimFastQC) {
            FASTQC(ch_reads_pre_align, "trimmed")
            MULTIQC(FASTQC.out.fastq_ch.collect(), "trimmed")
        }

    } else {
        ch_reads_pre_align = ch_readsRaw
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

}
