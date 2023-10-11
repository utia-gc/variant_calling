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
        tuple val(metadata), path(reads1), path(reads2)

    output:
        path('*.html'), emit: html
        path('*.zip'),  emit: zip

    script:
        if(metadata.readType == 'single') {
            """
            fastqc \
                --quiet \
                --threads ${task.cpus} \
                ${reads1}
            """
        } else if(metadata.readType == 'paired') {
            """
            fastqc \
                --quiet \
                --threads ${task.cpus} \
                ${reads1} \
                ${reads2}
            """
        }
}
