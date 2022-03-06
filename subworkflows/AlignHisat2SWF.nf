include { Hisat2Build } from '../modules/Hisat2Build.nf'

workflow AlignHisat2SWF {
    take:

    main:
        if (params.hisat2) {
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
}