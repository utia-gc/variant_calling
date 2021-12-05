process fastQC {
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(name), file(fastq_file)

    output:
    publishDir "${params.baseDirReport}/readsQC", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC", mode: 'copy', pattern: '*.zip'
    path "*", emit: ch_fastQC

    script:
    """
    fastqc ${fastq_file}
    """
}