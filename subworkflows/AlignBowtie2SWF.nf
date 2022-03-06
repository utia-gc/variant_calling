include { Bowtie2Build } from '../modules/Bowtie2Build.nf'
include { AlignBowtie2 } from '../modules/AlignBowtie2.nf'

workflow AlignBowtie2SWF {
    take:
        reads
    
    main:
        if (params.bowtie2) {
            bowtie2Indexes = Channel
                .fromPath("${params.bowtie2}*", checkIfExists: true)
                .collect()
        } else {
            Bowtie2Build(
                params.fasta
            )
            .collect()
            .set { bowtie2Indexes }
        }

        AlignBowtie2(reads, bowtie2Indexes)

    emit:
        sam = AlignBowtie2.out.sam
}
