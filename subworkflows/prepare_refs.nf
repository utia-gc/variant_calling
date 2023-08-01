include { GUNZIP as GUNZIP_FASTA       } from "../modules/gunzip.nf"
include { GUNZIP as GUNZIP_ANNOTATIONS } from "../modules/gunzip.nf"


workflow PREPARE_REFS {
    take:
        fasta
        annotations

    main:
        if(fasta.toString().endsWith(".gz")) {
            GUNZIP_FASTA(
                fasta
            )
            ch_fasta = GUNZIP_FASTA.out.gunzip
        } else {
            ch_fasta = Channel.value(file(fasta))
        }

        if(annotations.toString().endsWith(".gz")) {
            GUNZIP_ANNOTATIONS(
                annotations
            )
            ch_annotations = GUNZIP_ANNOTATIONS.out.gunzip
        } else {
            ch_annotations = Channel.value(file(annotations))
        }

    emit:
        fasta       = ch_fasta
        annotations = ch_annotations
}
