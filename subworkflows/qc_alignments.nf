include { multiqc as multiqc_alignments } from "../modules/multiqc.nf"
include { samtools_idxstats             } from '../modules/samtools_idxstats.nf'
include { samtools_stats                } from '../modules/samtools_stats.nf'


workflow QC_Alignments {
    take:
        alignments
        // map_log

    main:
        samtools_idxstats(alignments)
        samtools_stats(alignments)

        ch_multiqc_alignments = Channel.empty()
            .concat(samtools_idxstats.out.idxstat)
            .concat(samtools_stats.out.samtools_stats)
            .concat(samtools_stats.out.samtools_IS)
            .concat(samtools_stats.out.samtools_COV)
            // .concat(map_log)
            .collect( sort: true )

        multiqc_alignments(
            ch_multiqc_alignments,
            file(params.multiqcConfig),
            'alignments'
        )

    emit:
        multiqc = ch_multiqc_alignments
}
