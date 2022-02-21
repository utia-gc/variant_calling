include { FastQC } from '../modules/RawFastQC.nf'
include { ReadsMultiQC } from '../modules/RawReadsMultiQC.nf'

workflow RawReadsQCSWF {
    take:
        rawReads
        runName

    main:
        FastQC(rawReads)
        ReadsMultiQC(FastQC.out.zip.collect { it[1] } , runName)

    emit:
        fqc_zip = FastQC.out.zip
}
