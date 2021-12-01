process multiQC {
    container 'ewels/multiqc:v1.11'

    input:
    path fastqc

    output:
    publishDir "${params.baseDirReport}/fastQC", mode: 'copy'
    path "*"

    script:
    """
    multiqc ${fastqc}
    """
}