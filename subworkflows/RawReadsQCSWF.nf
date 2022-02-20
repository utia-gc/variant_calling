include { FastQC } from '../modules/FastQC.nf'

workflow RawReadsQCSWF {
    take:
        rawReads

    main:
        FastQC(rawReads)

    emit:
        fqc_zip = FastQC.out.zip
}