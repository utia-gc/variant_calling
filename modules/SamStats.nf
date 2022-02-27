/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-27
Purpose: Nextflow module
*/

process SamStats {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/samtools:1.15--h1170115_1'

    publishDir "${params.baseDirData}/align/stats", mode: 'copy', pattern: '*.txt'

    input:
        tuple val(metadata), path(bam), path(bai)

    output:
        tuple val(metadata), path("*_sST*"), path("*_sIX*"), emit: samStats

    script:
        """
        samtools stats ${bam} > ${metadata.sampleName}__sbl_sSR_sST.txt
        samtools idxstats ${bam} > ${metadata.sampleName}__sbl_sSR_sIX.txt
        """
}
