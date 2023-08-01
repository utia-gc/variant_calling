include { gunzip as gunzip_fasta      } from "../modules/gunzip.nf"
include { gunzip as gunzip_annotations } from "../modules/gunzip.nf"


workflow PREPARE_REFS {
    take:
        fasta
        annotations

    main:
        if(fasta.toString().endsWith(".gz")) {
            gunzip_fasta(
                fasta
            )
            ch_fasta = gunzip_fasta.out.gunzip
        } else {
            ch_fasta = Channel.value(file(fasta))
        }

        if(annotations.toString().endsWith(".gz")) {
            gunzip_annotations(
                annotations
            )
            ch_annotations = gunzip_annotations.out.gunzip
        } else {
            ch_annotations = Channel.value(file(annotations))
        }

    emit:
        fasta       = ch_fasta
        annotations = ch_annotations
}
