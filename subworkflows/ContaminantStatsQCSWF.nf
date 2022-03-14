include { SamtoolsFlagstat   } from '../modules/SamtoolsFlagstat.nf'
include { ContaminantStatsQC } from '../modules/ContaminantStatsQC.nf'

workflow ContaminantStatsQCSWF {
    take:
        sam
        runName

    main:
        SamtoolsFlagstat(
            sam
        )
        ContaminantStatsQC(
            SamtoolsFlagstat.out.sFS.collect(),
            runName,
            SamtoolsFlagstat.out.tools.first()
        )
    
    emit:
        samtoolsFlagstat = SamtoolsFlagstat.out.sFS
}