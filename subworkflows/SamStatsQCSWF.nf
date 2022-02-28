include { SamStats } from '../modules/SamStats.nf'
include { SamStatsMultiQC } from '../modules/SamStatsMultiQC.nf'

workflow SamStatsQCSWF {
    take:
        bamBai
        runName
    
    main:
        SamStats(bamBai)
        SamStatsMultiQC(SamStats.out.sST.collect { it[1] } , SamStats.out.sIX.collect { it[1] } , runName)
}
