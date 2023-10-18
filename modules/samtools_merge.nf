process samtools_merge {
    tag "${metadata.sampleName}"

    label 'samtools'

    input:
        tuple val(metadata), path(bams), path(bais)

    output:
        tuple val(metadata), path('*_merged.bam'), emit: bamMerged

    script:
        """
        samtools merge \
            -o ${metadata.sampleName}_merged.bam \
            ${bams}
        """
}