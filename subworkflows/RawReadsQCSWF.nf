include { FastQC } from '../modules/FastQC.nf'
include { ReadsMultiQC } from '../modules/ReadsMultiQC.nf'

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
