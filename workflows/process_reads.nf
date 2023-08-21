include { Cat_Reads  } from "../subworkflows/cat_reads.nf"
include { Trim_Reads } from "../subworkflows/trim_reads.nf"

/**
 * Workflow to process reads.
 * Reads should be trimmed and concatenated here, if necessary.
 */

workflow PROCESS_READS {
    take:
        reads_raw

    main:
        Cat_Reads(reads_raw)

        if(!params.skipTrimReads) {
            Trim_Reads(
                Cat_Reads.out.reads_catted,
                params.tools.trim
            )
            ch_reads_post_trim = Trim_Reads.out.reads_trim
            ch_trim_log        = Trim_Reads.out.trim_log
        } else {
            ch_reads_post_trim = Cat_Reads.out.reads_catted
            ch_trim_log        = Channel.empty()
        }

    emit:
        reads_pre_align = ch_reads_post_trim
        trim_log        = ch_trim_log .ifEmpty('EMPTY')
}
