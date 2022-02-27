include { AlignBowtie2 } from '../modules/AlignBowtie2.nf'
include { CompressSortSam } from '../modules/CompressSortSam.nf'

workflow AlignBowtie2SWF {
    take:
        reads
        runName
        bt2IndexBase
    
    main:
        bt2Indexes = Channel
            .fromPath("${bt2IndexBase}*", checkIfExists: true)
            .collect()

        AlignBowtie2(reads, bt2Indexes)

        CompressSortSam(AlignBowtie2.out.sam)

    emit:
        sam = AlignBowtie2.out.sam
        bam = CompressSortSam.out.bam
}
