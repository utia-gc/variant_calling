include { gunzip as gunzip_genome      } from "../modules/gunzip.nf"
include { gunzip as gunzip_annotations } from "../modules/gunzip.nf"


workflow Prepare_Refs {
    take:
        genome
        annotations

    main:
        if(genome.toString().endsWith(".gz")) {
            gunzip_genome(
                genome
            )
            ch_genome = gunzip_genome.out.gunzip
        } else {
            ch_genome = Channel.value(file(genome))
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
        genome      = ch_genome
        annotations = ch_annotations
}
