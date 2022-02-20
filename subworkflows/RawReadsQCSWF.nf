include { FastQC } from '../modules/FastQC.nf'
include { ReadsMultiQC } from '../modules/ReadsMultiQC.nf'

workflow RawReadsQCSWF {
    take:
        rawReads

    main:
        FastQC(rawReads)
        ReadsMultiQC(FastQC.out.zip.collect { it[1] } )

    emit:
        fqc_zip = FastQC.out.zip
}
