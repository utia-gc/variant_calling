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
        path '*_samtools-IS.txt',    emit: samtools_IS

    script:
        """
        samtools stats \
            ${bam} \
            > ${metadata.sampleName}_samtools-stats.txt

        # extract insert sizes
        grep ^IS ${metadata.sampleName}_samtools-stats.txt \
            | cut -f 2,3 \
            > ${metadata.sampleName}_samtools-IS.txt
        """
}
