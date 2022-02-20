process ReadsMultiQC {
    container 'ewels/multiqc:v1.11'

    publishDir "${params.baseDirReport}/readsQC", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC", mode: 'copy', pattern: '*multiqc_data*'

    input:
        file(fastqc)

    output:
        path "*"

    script:
        """
        multiqc ${fastqc}
        """
}
