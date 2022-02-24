include { FastQC } from '../modules/TrimFastQC.nf'
include { ReadsMultiQC } from '../modules/TrimReadsMultiQC.nf'

workflow TrimReadsQCSWF {
    take:
        trimReads
        runName

    main:
        FastQC(trimReads)
        ReadsMultiQC(FastQC.out.zip.collect { it[1] } , runName)

    emit:
        fqc_zip = FastQC.out.zip
}
