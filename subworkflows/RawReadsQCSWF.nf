include { FastQC } from '../modules/RawFastQC.nf'
include { ReadsMultiQC } from '../modules/RawReadsMultiQC.nf'

workflow RawReadsQCSWF {
    take:
        rawReads
        runName

    main:
        FastQC(rawReads)
        ReadsMultiQC(FastQC.out.zip.collect(), runName, FastQC.out.tools.first())

    emit:
        fqc_zip = FastQC.out.zip
}
