include { fastqc as fastqc_raw        } from "../modules/fastqc.nf"
include { fastqc as fastqc_prealign   } from "../modules/fastqc.nf"
include { multiqc as multiqc_raw      } from "../modules/multiqc.nf"
include { multiqc as multiqc_prealign } from "../modules/multiqc.nf"
include { multiqc as multiqc_full     } from "../modules/multiqc.nf"

workflow CHECK_QUALITY {
    take:
        reads_raw
        reads_prealign
        trim_log

    main:
        if(!params.skipRawFastQC) {
            fastqc_raw(reads_raw)
            ch_multiqc_reads_raw = Channel.empty()
                .concat(fastqc_raw.out.zip)
                .collect( sort: true )
        } else {
            ch_multiqc_reads_raw = Channel.empty()
        }
        multiqc_raw(
            ch_multiqc_reads_raw,
            "raw"
        )

        if(!params.skipTrimFastQC) {
            fastqc_prealign(reads_prealign)
            ch_multiqc_reads_prealign = Channel.empty()
                .concat(fastqc_prealign.out.zip)
                .concat(trim_log)
                .collect( sort: true )
        } else {
            ch_multiqc_reads_prealign = Channel.empty()
        }
        multiqc_prealign(
            ch_multiqc_reads_prealign,
            "prealign"
        )

        ch_multiqc_full = Channel.empty()
            .concat(ch_multiqc_reads_raw)
            .concat(ch_multiqc_reads_prealign)
            .collect( sort: true )
        multiqc_full(
            ch_multiqc_full,
            "full"
        )
}
