include { Bowtie2Build } from '../modules/Bowtie2Build.nf'
include { Bowtie2Align } from '../modules/Bowtie2Align.nf'

workflow AlignBowtie2SWF {
    take:
        reads
        index
    
    main:
        if (index) {
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

        Bowtie2Align(reads, bowtie2Indexes)

    emit:
        sam = Bowtie2Align.out.sam
}
