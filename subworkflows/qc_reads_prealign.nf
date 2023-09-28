include { fastqc           as fastqc_prealign  } from '../modules/fastqc.nf'
include { multiqc          as multiqc_prealign } from '../modules/multiqc.nf'
include { sequencing_depth                     } from '../modules/sequencing_depth.nf'


workflow QC_Reads_Prealign {
    take:
        reads_prealign
        trim_log
        genome_index

    main:
        fastqc_prealign(reads_prealign)

        sequencing_depth(
            reads_prealign,
            genome_index
        )

        ch_multiqc_reads_prealign = Channel.empty()
            .concat(fastqc_prealign.out.zip)
            .concat(sequencing_depth.out.depth)
            .concat(trim_log)
            .collect()

        multiqc_prealign(
            ch_multiqc_reads_prealign,
            file(params.multiqcConfig),
            "reads_prealign"
        )

    emit:
        multiqc = ch_multiqc_reads_prealign
}
