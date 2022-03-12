include { Hisat2Build } from '../modules/Hisat2Build.nf'
include { Hisat2Align } from '../modules/Hisat2Align.nf'

workflow AlignHisat2SWF {
    take:
        reads
        index

    main:
        // set or build hisat2 index
        if (index) {
            hisat2Indexes = Channel
                .fromPath("${params.hisat2}*", checkIfExists: true)
                .collect()
        } else {
            Hisat2Build(
                params.fasta
            )
            .collect()
            .set { hisat2Indexes }
        }

        Hisat2Align(
            reads,
            hisat2Indexes
        )

    emit:
        sam = Hisat2Align.out.sam
}