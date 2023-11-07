include { gatk_HaplotypeCaller } from '../modules/gatk_HaplotypeCaller.nf'


workflow GENOTYPE {
    take:
        alignments
        genome

    main:
        gatk_HaplotypeCaller(
            alignments,
            genome
        )
}
