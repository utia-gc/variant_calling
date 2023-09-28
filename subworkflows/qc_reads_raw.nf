include { fastqc           as fastqc_raw  } from "../modules/fastqc.nf"
include { multiqc          as multiqc_raw } from "../modules/multiqc.nf"
include { sequencing_depth                } from '../modules/sequencing_depth.nf'

workflow QC_Reads_Raw {
    take:
        reads_raw
        genome_index

    main:
        fastqc_raw(reads_raw)

        sequencing_depth(
            reads_raw,
            genome_index
        )

        ch_multiqc_reads_raw = Channel.empty()
            .concat(fastqc_raw.out.zip)
            .concat(sequencing_depth.out.depth)
            .collect()

        multiqc_raw(
            ch_multiqc_reads_raw,
            file(params.multiqcConfig),
            "reads_raw"
        )

    emit:
        multiqc = ch_multiqc_reads_raw
}
