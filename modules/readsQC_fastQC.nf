process fastQC {
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(name), file(fastq_file)

    output:
    publishDir "${params.baseDirReport}/fastQC", mode: 'copy'
    path "*", emit: ch_fastQC

    script:
    """
    fastqc ${fastq_file}
    """
}