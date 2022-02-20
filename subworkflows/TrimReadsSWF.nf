include { Fastp } from '../modules/Fastp.nf'

workflow TrimReadsSWF {
    take:
        rawReads

    main:
        Fastp(rawReads)
}