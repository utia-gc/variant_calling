include { cutadapt } from "../modules/cutadapt.nf"

/**
 * Workflow to process reads.
 * Reads should be trimmed and concatenated here, if necessary.
 */

workflow PROCESS_READS {
    take:
        reads_raw

    main:
        if(!params.skipTrimReads) {
            cutadapt(
                reads_raw,
                params.r1_adapter,
                params.r2_adapter,
                params.minimum_length
            )
            ch_reads_trim = cutadapt.out.reads
            ch_trim_log   = cutadapt.out.log
        } else {
            ch_reads_trim = Channel.empty()
            ch_trim_log   = Channel.empty()
        }

    emit:
        reads_pre_align = ch_reads_trim
        trim_log        = ch_trim_log
}
