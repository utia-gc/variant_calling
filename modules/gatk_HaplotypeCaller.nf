/**
 * Process to run GATK HaplotypeCaller
 *
 * Call germline SNPs and indels via local re-assembly of haplotypes.
 * @see https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller
 *
 * @input alignments the alignments channel of format [metadata, bam, bam.bai]
 * @input reference the reference sequence file channel of format [fasta, fai, dict]
 * @emit gvcf the output GVCF file channel of format [metadata, GVCF]. The GVCF is gzipped.
 */
process gatk_HaplotypeCaller {
    tag "${metadata.sampleName}"

    label 'gatk'

    label 'big_cpu'
    label 'big_mem'
    label 'max_time'

    input:
        tuple val(metadata), path(bam), path(bai)
        tuple path(fasta), path(fai), path(dict)

    output:
        tuple val(metadata), path('*.g.vcf.gz'), emit: gvcf

    script:
        String args = new Args(task.ext).buildArgsString()

        """
        gatk HaplotypeCaller \\
            --input ${bam} \\
            --output ${metadata.sampleName}.g.vcf.gz \\
            --reference ${fasta} \\
            --emit-ref-confidence GVCF \\
            ${args}
        """
}
