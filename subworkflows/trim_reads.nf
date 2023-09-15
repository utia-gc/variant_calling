include { cutadapt } from "../modules/cutadapt.nf"
include { fastp    } from "../modules/fastp.nf"

/**
 * Workflow to process reads.
 * Reads should be trimmed and concatenated here, if necessary.
 */

workflow Trim_Reads {
    take:
        reads_raw
        trim_tool

    main:
        switch( trim_tool.toUpperCase() ) {
            case "CUTADAPT":
                cutadapt(
                    reads_raw,
                    params.r1_adapter,
                    params.r2_adapter,
                    params.minimum_length
                )
                ch_reads_trim = cutadapt.out.reads
                ch_trim_log   = cutadapt.out.log
                break

            case "FASTP":
                fastp(reads_raw)
                ch_reads_trim = fastp.out.reads
                ch_trim_log   = fastp.out.log
                break
        }

    emit:
        reads_trim = ch_reads_trim
        trim_log   = ch_trim_log
}
