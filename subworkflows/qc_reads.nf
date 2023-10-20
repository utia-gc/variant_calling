include { fastqc  as fastqc_prealign  } from '../modules/fastqc.nf'
include { fastqc  as fastqc_raw       } from "../modules/fastqc.nf"
include { multiqc as multiqc_reads    } from '../modules/multiqc.nf'

workflow QC_Reads {
    take:
        reads_raw
        reads_prealign
        trim_log

    main:
        fastqc_raw(reads_raw)
        fastqc_prealign(reads_prealign)

        ch_multiqc_reads = Channel.empty()
            .concat(fastqc_raw.out.zip)
            .concat(fastqc_prealign.out.zip)
            .concat(trim_log)
            .collect()

        multiqc_reads(
            ch_multiqc_reads,
            file("${projectDir}/assets/multiqc_config.yaml"),
            "reads"
        )

    emit:
        multiqc = ch_multiqc_reads
}
