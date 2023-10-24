process gatk_MarkDuplicates {
    tag "${metadata.sampleName}"

    label 'gatk'

    publishDir(
        path:    "${params.publishDirData}/alignments",
        mode:    "${params.publishMode}",
        pattern: "${metadata.sampleName}.bam{,.bai}"
    )

    input:
        tuple val(metadata), path(bam)

    output:
        tuple val(metadata), path('*.bam'), path('*.bam.bai'), emit: bamMarkDupIndexed

    script:
        """
        gatk MarkDuplicates \
            --INPUT ${bam} \
            --METRICS_FILE ${metadata.sampleName}_MarkDuplicates-metrics.txt \
            --OUTPUT ${metadata.sampleName}.bam \
            --CREATE_INDEX \
            --VALIDATION_STRINGENCY SILENT

        mv ${metadata.sampleName}.bai ${metadata.sampleName}.bam.bai
        """
}
