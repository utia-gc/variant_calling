process gatk_MergeSameFiles {
    tag "${metadata.sampleName}"

    label 'gatk'

    input:
        tuple val(metadata), path(bams), path(bais)

    output:
        tuple val(metadata), path('*_merged.bam'), emit: bamMerged

    script:
        // create the string of INPUT arguments
        // prepend '--INPUT' to each BAM
        // join into a space-delimted string
        def inputs = bams.collect { bam ->
            "--INPUT ${bam}" 
        }.join(' ')

        """
        gatk MergeSamFiles \
            ${inputs} \
            --OUTPUT ${metadata.sampleName}_merged.bam \
            --VALIDATION_STRINGENCY SILENT
        """
}
