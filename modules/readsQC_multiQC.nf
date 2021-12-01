process multiQC {
    container 'ewels/multiqc:v1.11'

    input:
    path fastqc

    output:
    publishDir "${params.baseDirReport}/readsQC", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC", mode: 'copy', pattern: '*multiqc_data*'
    path "*"

    script:
    """
    multiqc ${fastqc}
    """
}