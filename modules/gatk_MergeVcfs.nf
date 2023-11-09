/**
 * Process to run GATK MergeVcfs
 *
 * Combines multiple variant files into a single variant file.
 * @see https://gatk.broadinstitute.org/hc/en-us/articles/13832753837467-MergeVcfs-Picard-
 *
 * @input vcfs the joint genotype files for each interval in gzipped VCF format.
 * @emit vcf the merged joint genotype file in gzipped VCF format and its index file of format [vcf.gz, vcf.gz.tbi].
 */
process gatk_MergeVcfs {
    label 'gatk'

    publishDir(
        path:    "${params.publishDirData}/variants",
        mode:    "${params.publishMode}",
        pattern: "*.vcf.gz{,.tbi}"
    )

    input:
        path vcfs

    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: vcf

    script:
        // create the string of INPUT arguments
        // prepend '--INPUT' to each VCF
        // join into a space-delimted string
        def inputs = vcfs.collect { vcf ->
            "--INPUT ${vcf}" 
        }.join(' ')

        String args = new Args(task.ext).buildArgsString()

        """
        mkdir \${PWD}/tmp

        gatk MergeVcfs \\
            ${inputs} \\
            --OUTPUT merged.vcf.gz \\
            ${args}

        gatk IndexFeatureFile \\
            --input merged.vcf.gz
        """
}
