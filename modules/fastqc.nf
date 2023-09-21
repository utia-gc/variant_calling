process fastqc {
    tag "${metadata.sampleName}"
    
    label 'fastqc'

    label 'lil_mem'

    publishDir(
        path:    "${params.publishDirReports}/fastqc",
        mode:    "${params.publishMode}",
        pattern: '*{.html,.zip}'
    )

    input:
        tuple val(metadata), path(reads)

    output:
        path('*.html'), emit: html
        path('*.zip'),  emit: zip

    script:
        """
        fastqc \
            --quiet \
            --threads ${task.cpus} \
            ${reads}
        """
}
