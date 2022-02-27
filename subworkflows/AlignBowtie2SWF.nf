include { AlignBowtie2    } from '../modules/AlignBowtie2.nf'
include { Samblaster      } from '../modules/Samblaster.nf'
include { CompressSortSam } from '../modules/CompressSortSam.nf'
include { IndexBam        } from '../modules/IndexBam.nf'
include { SamStats        } from '../modules/SamStats.nf'

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

        Samblaster(AlignBowtie2.out.sam)

        CompressSortSam(Samblaster.out.sam)

        IndexBam(CompressSortSam.out.bam)

        SamStats(IndexBam.out.bamBai)

    emit:
        bam = CompressSortSam.out.bam
        bamBai = IndexBam.out.bamBai
}
