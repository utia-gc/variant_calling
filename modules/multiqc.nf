process multiqc {
    tag "${fileName}"
    
    label 'multiqc'

    label 'med_mem'

    publishDir(
        path:    "${params.publishDirReports}/multiqc/${fileName}",
        mode:    "${params.publishMode}"
    )

    input:
        path('*')
        path config
        val fileName

    output:
        path("*html"),              emit: report
        path("${fileName}_data/*"), emit: data

    script:
        """
        multiqc \
            --filename ${fileName} \
            --config ${config} \
            .
        """
}
