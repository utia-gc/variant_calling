include { AlignBowtie2 } from '../modules/AlignBowtie2.nf'

workflow AlignBowtie2SWF {
    take:
        reads
        runName
        bt2Indexes
    
    main:
        AlignBowtie2(reads, runName, bt2Indexes)

    emit:
        sam = AlignBowtie2.out.sam
}