include { fastqc  as fastqc_raw  } from "../modules/fastqc.nf"
include { multiqc as multiqc_raw } from "../modules/multiqc.nf"

workflow QC_Reads_Raw {
    take:
        reads_raw

    main:
        fastqc_raw(reads_raw)

        ch_multiqc_reads_raw = Channel.empty()
            .concat(fastqc_raw.out.zip)
            .collect()

        multiqc_raw(
            ch_multiqc_reads_raw,
            file("${projectDir}/assets/multiqc_config.yaml"),
            "reads_raw"
        )

    emit:
        multiqc = ch_multiqc_reads_raw
}
