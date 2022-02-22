process ReadsMultiQC {
    container 'ewels/multiqc:v1.11'

    publishDir "${params.baseDirReport}/readsQC/trim", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC/trim", mode: 'copy', pattern: '*multiqc_data*'

    input:
        file(fastqc)
        val(runName)

    output:
        path "*"

    script:
        """
        multiqc -n ${runName}__fsp_fqc_mqc ${fastqc}
        """
}
