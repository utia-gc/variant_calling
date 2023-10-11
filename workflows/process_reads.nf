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
        log.warn "Concatenate reads is deprecated and is no longer the default behavior."
        log.warn "If you wish to concatenate reads, you must specifically set `params.concatenateReads = true`"
        if(params.concatenateReads) {
            Cat_Reads(reads_raw)
            ch_reads_post_cat = Cat_Reads.out.reads_catted
        } else {
            ch_reads_post_cat = reads_raw
        }

        if(!params.skipTrimReads) {
            Trim_Reads(
                ch_reads_post_cat,
                params.tools.trim
            )
            ch_reads_post_trim = Trim_Reads.out.reads_trim
            ch_trim_log        = Trim_Reads.out.trim_log
        } else {
            ch_reads_post_trim = ch_reads_post_cat
            ch_trim_log        = Channel.empty()
        }

    emit:
        reads_pre_align = ch_reads_post_trim
        trim_log        = ch_trim_log
}
