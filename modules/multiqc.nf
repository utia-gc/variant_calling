process multiqc {
    tag "${fileName}"
    
    label 'multiqc'
    label 'lil_mem'

    publishDir(
        path:    "${params.publishDirReports}/multiqc/${fileName}",
        mode:    "${params.publishMode}"
    )

    input:
        path('*')
        val fileName

    output:
        path("*html"),              emit: report
        path("${fileName}_data/*"), emit: data

    script:
        """
        multiqc \
            --filename ${fileName} \
            .
        """
}
