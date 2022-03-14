include { Preseq           } from '../modules/Preseq.nf'
include { PreseqRealCounts } from '../modules/PreseqRealCounts.nf'

workflow PreseqSWF {
    take:
        bam
    
    main:
        Preseq(
            bam
        )
        PreseqRealCounts(
            bam,
            Preseq.out.psL
        )
    
    emit:
        psL = Preseq.out.psL
        psRealCounts = PreseqRealCounts.out.preseqRealCounts
}