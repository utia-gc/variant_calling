include { fastqc as fastqc_prealign   } from "../modules/fastqc.nf"
include { multiqc as multiqc_prealign } from "../modules/multiqc.nf"

workflow QC_Reads_Prealign {
    take:
        reads_prealign
        trim_log

    main:
        fastqc_prealign(reads_prealign)

        ch_multiqc_reads_prealign = Channel.empty()
            .concat(fastqc_prealign.out.zip)
            .concat(trim_log)
            .collect( sort: true )

        multiqc_prealign(
            ch_multiqc_reads_prealign,
            file(params.multiqcConfig),
            "reads_prealign"
        )

    emit:
        multiqc = ch_multiqc_reads_prealign
}
