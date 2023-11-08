/**
 * Index a GVCF file with GATK IndexFeatureFile
 *
 * Call germline SNPs and indels via local re-assembly of haplotypes.
 * @see https://gatk.broadinstitute.org/hc/en-us/articles/13832767900443-IndexFeatureFile
 *
 * @input gvcf the GVCF channel of format [metadata, g.vcf.gz] where GVCF is gzipped.
 * @emit gvcfTbi the output GVCF file channel of format [metadata, g.vcf.gz, g.vcf.gz.tbi].
 */
process gatk_IndexGVCF {
    tag "${metadata.sampleName}"

    label 'gatk'

    input:
        tuple val(metadata), path(gvcf)

    output:
        tuple val(metadata), path(gvcf), path('*.g.vcf.gz.tbi'), emit: gvcfTbi

    script:
        String args = new Args(task.ext).buildArgsString()

        """
        gatk IndexFeatureFile \\
            --input ${metadata.sampleName}.g.vcf.gz
        """
}
