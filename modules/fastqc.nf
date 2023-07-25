process FASTQC {
    tag "${metadata.sampleName}"
    
    label 'fastqc'
    label 'lil_mem'

    publishDir(path: "${publish_dir}/fastqc", mode: "copy")

    input:
        tuple val(metadata), path(reads)

    output:
        path('*.html'), emit: html
        path('*.zip'),  emit: zip

    script:
        """
        fastqc \
            -q \
            ${reads}
        """
}
