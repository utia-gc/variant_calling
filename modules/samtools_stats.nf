process samtools_stats {
    tag "${metadata.sampleName}"

    label 'samtools'

    publishDir(
        path:    "${params.publishDirReports}/alignments",
        mode:    "${params.publishMode}",
        pattern: '*_samtools-stats.txt'
    )

    input:
        tuple val(metadata), path(bam), path(bai)

    output:
        path '*_samtools-stats.txt', emit: samtools_stats

    script:
        """
        samtools stats \
            ${bam} \
            > ${metadata.sampleName}_samtools-stats.txt
        """
}
