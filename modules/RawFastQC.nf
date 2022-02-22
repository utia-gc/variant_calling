process FastQC {
    tag "${metadata.sampleName}"

    container 'biocontainers/fastqc:v0.11.9_cv8'

    publishDir "${params.baseDirReport}/readsQC/raw/fastQC", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/readsQC/raw/fastQC",   mode: 'copy', pattern: '*.zip'

    input:
        tuple val(metadata), file(reads)

    output:
        tuple val(metadata), path('*.html'), emit: html
        tuple val(metadata), path('*.zip'), emit: zip 

    script:
        if (metadata.readType == 'single') {
            read1 = "${metadata.sampleName}_R1.fastq.gz"

            """
            [ ! -f ${read1} ] && ln -s ${reads} ${read1}

            fastqc ${read1}
            """
        } else {
            read1 = "${metadata.sampleName}_R1.fastq.gz"
            read2 = "${metadata.sampleName}_R2.fastq.gz"

            """
            [ ! -f ${read1} ] && ln -s ${reads[0]} ${read1}
            [ ! -f ${read2} ] && ln -s ${reads[1]} ${read2}

            fastqc ${read1} ${read2}
            """
        }
}
