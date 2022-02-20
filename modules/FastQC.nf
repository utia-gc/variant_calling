process FastQC {
    tag "${metadata.sampleName}"

    container 'biocontainers/fastqc:v0.11.9_cv8'

    publishDir "${params.baseDirReport}/readsQC/fastQC", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC/fastQC",   mode: 'copy', pattern: '*.zip'

    input:
        tuple val(metadata), file(reads)

    output:
        tuple val(metadata), path('*.html'), emit: html
        tuple val(metadata), path('*.zip'), emit: zip 

    script:
        """
        fastqc ${reads}
        """
}
