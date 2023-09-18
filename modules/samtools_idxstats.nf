process samtools_idxstats {
    tag "${metadata.sampleName}"

    label 'samtools'

    publishDir(
        path:    "${params.publishDirReports}/alignments",
        mode:    "${params.publishMode}",
        pattern: '*_idxstat.txt'
    )

    input:
        tuple val(metadata), path(bam), path(bai)

    output:
        path '*_idxstat.txt', emit: idxstat

    script:
        """
        samtools idxstats \
            ${bam} \
            > ${metadata.sampleName}_idxstat.txt
        """
}
