/**
 * Process to run GATK GenotypeGVCFs
 *
 * Perform joint genotyping on one or more samples pre-called with HaplotypeCaller and stored in GenomicsDB.
 * @see https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs
 *
 * @input genomicsDB the GenomicsDB of merged GVCF files and interval used of foramt [interval, GenomicsDB].
 * @input reference the reference sequence file channel of format [fasta, fai, dict].
 * @emit vcf the joint genotype file in gzipped VCF format.
 */
process gatk_GenotypeGVCFs {
    label 'gatk'

    input:
        tuple val(interval), path(genomicsDB)
        tuple path(fasta), path(fai), path(dict)

    output:
        path '*.vcf.gz', emit: vcf

    script:
        String args = new Args(task.ext).buildArgsString()

        String javaOptions = ""

        """
        mkdir \${PWD}/tmp

        gatk --java-options "${javaOptions}" GenotypeGVCFs \\
            --variant gendb://${genomicsDB} \\
            --reference ${fasta} \\
            --output ${interval}.vcf.gz \\
            --intervals ${interval} \\
            ${args}
        """
}
