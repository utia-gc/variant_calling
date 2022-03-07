include { Samblaster      } from '../modules/Samblaster.nf'
include { CompressSortSam } from '../modules/CompressSortSam.nf'
include { IndexBam        } from '../modules/IndexBam.nf'

workflow PreprocessSamSWF {
    take:
        sam

    main:
        // mark duplicates with samblaster
        Samblaster(
            sam
        )

        // compress sam and sort bam by coordinate
        CompressSortSam(
            Samblaster.out.sam
        )

        // index bam
        IndexBam(
            CompressSortSam.out.bam
        )

    emit:
        bam    = CompressSortSam.out.bam
        bamBai = IndexBam.out.bamBai
}