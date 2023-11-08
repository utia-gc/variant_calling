include { gatk_HaplotypeCaller } from '../modules/gatk_HaplotypeCaller.nf'
include { gatk_IndexGVCF       } from '../modules/gatk_IndexGVCF.nf'


workflow GENOTYPE {
    take:
        alignments
        genome

    main:
        gatk_HaplotypeCaller(
            alignments,
            genome
        )
          | gatk_IndexGVCF
}
