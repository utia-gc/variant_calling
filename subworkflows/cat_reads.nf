include { cat_fastq   } from '../modules/cat_fastq.nf'
include { Group_Reads } from '../subworkflows/group_reads.nf'

workflow Cat_Reads {
    take:
        reads

    main:
        reads \
            | Group_Reads \
            | cat_fastq

    emit:
        reads_catted = cat_fastq.out.reads
}
