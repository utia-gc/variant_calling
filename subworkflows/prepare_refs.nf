include { gatk_CreateSequenceDictionary        } from '../modules/gatk_CreateSequenceDictionary.nf'
include { gunzip         as gunzip_genome      } from "../modules/gunzip.nf"
include { gunzip         as gunzip_annotations } from "../modules/gunzip.nf"
include { samtools_faidx                       } from '../modules/samtools_faidx.nf'


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

        samtools_faidx(ch_genome)
        ch_genome_index = samtools_faidx.out.fai

        gatk_CreateSequenceDictionary(ch_genome)
        ch_genome_dict = gatk_CreateSequenceDictionary.out.dict

        if(annotations.toString().endsWith(".gz")) {
            gunzip_annotations(
                annotations
            )
            ch_annotations = gunzip_annotations.out.gunzip
        } else {
            ch_annotations = Channel.value(file(annotations))
        }

    emit:
        genome       = ch_genome
        genome_index = ch_genome_index
        genome_dict  = ch_genome_dict
        annotations  = ch_annotations
}
