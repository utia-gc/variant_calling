include { gatk_GenomicsDBImport } from '../modules/gatk_GenomicsDBImport.nf'
include { gatk_GenotypeGVCFs    } from '../modules/gatk_GenotypeGVCFs.nf'
include { gatk_HaplotypeCaller  } from '../modules/gatk_HaplotypeCaller.nf'
include { gatk_IndexGVCF        } from '../modules/gatk_IndexGVCF.nf'


workflow GENOTYPE {
    take:
        alignments
        genome

    main:
        // call haplotypes and index GVCFs
        gatk_HaplotypeCaller(
            alignments,
            genome
        )
          | gatk_IndexGVCF

        // collect channels of GVCFs and indexes
        gatk_IndexGVCF.out.gvcfTbi
            .multiMap { metadata, gvcf, gvcfTbi ->
                gvcf:    gvcf
                gvcfTbi: gvcfTbi
            }
            .set { ch_gvcf_index }
        ch_gvcf_index.gvcf
            .collect()
            .set { ch_gvcfs }
        ch_gvcf_index.gvcfTbi
            .collect()
            .set { ch_gvcfTbis }

        // extract genome intervals from fasta file
        genome
            .map { fasta, fai, dict ->
                fasta
            }
            .splitFasta( record: [id: true] )
            .filter { record ->
                record.id =~ params.intervalFilter
            }
            .map { record ->
                record.id
            }
            .set { ch_intervals }

        gatk_GenomicsDBImport(
            ch_gvcfs,
            ch_gvcfTbis,
            genome,
            ch_intervals
        )

        gatk_GenotypeGVCFs(
            gatk_GenomicsDBImport.out.genomicsDB,
            genome
        )
}
