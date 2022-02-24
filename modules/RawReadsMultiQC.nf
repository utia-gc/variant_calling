process ReadsMultiQC {
    container 'ewels/multiqc:v1.11'

    publishDir "${params.baseDirReport}/readsQC/raw", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC/raw", mode: 'copy', pattern: '*multiqc_data*'

    input:
        file(fastqc)
        val(runName)

    output:
        path "*"

    script:
        """
        multiqc -n ${runName}__fqc_mqc ${fastqc}
        """
}
