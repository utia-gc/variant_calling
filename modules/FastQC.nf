process FastQC {
    tag "${metadata.sampleName}"

    container 'biocontainers/fastqc:v0.11.9_cv8'

    publishDir "${params.baseDirReport}/readsQC/fastQC", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC/fastQC",   mode: 'copy', pattern: '*.zip'

    input:
        tuple val(metadata), file(reads)

    output:
        tuple val(metadata), path('*.html'), emit: html
        tuple val(metadata), path('*.zip'), emit: zip 

    script:
        if (metadata.readType == 'single') {
            """
            [ ! -f ${metadata.sampleName}_R1.fastq.gz ] && ln -s ${reads} ${metadata.sampleName}_R1.fastq.gz

            fastqc ${metadata.sampleName}_R1.fastq.gz
            """
        } else {
            """
            [ ! -f ${metadata.sampleName}_R1.fastq.gz ] && ln -s ${reads[0]} ${metadata.sampleName}_R1.fastq.gz
            [ ! -f ${metadata.sampleName}_R2.fastq.gz ] && ln -s ${reads[1]} ${metadata.sampleName}_R2.fastq.gz

            fastqc ${metadata.sampleName}_R1.fastq.gz ${metadata.sampleName}_R2.fastq.gz
            """
        }
}
