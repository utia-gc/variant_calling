include { SamStats } from '../modules/SamStats.nf'

workflow SamStatsQCSWF {
    take:
        bamBai
        runName
    
    main:
        SamStats(bamBai)
}
