include { AlignBowtie2 } from '../modules/AlignBowtie2.nf'

workflow AlignBowtie2SWF {
    take:
        reads
        runName
        bt2IndexBase
    
    main:
        bt2Indexes = Channel
            .fromPath("${bt2IndexBase}*", checkIfExists: true)
            .collect()

        AlignBowtie2(reads, runName, bt2Indexes)

    emit:
        sam = AlignBowtie2.out.sam
}
