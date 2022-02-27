/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-27
Purpose: Nextflow module
*/

process CompressSortSam {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/samtools:1.15--h1170115_1'

    publishDir "${params.baseDirData}/align", mode: 'copy', pattern: '*.bam'

    input:
        tuple val(metadata), file(sam)

    output:
        tuple val(metadata), file('*.bam'), emit: bam

    script:
        """
        samtools view -bh ${sam} | \
        samtools sort -o ${metadata.sampleName}__sSR.bam -
        """
}
