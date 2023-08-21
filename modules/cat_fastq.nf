process cat_fastq {
    tag "${metadata.sampleName}"

    label 'base'

    input:
        tuple val(metadata), path(reads1), path(reads2)

    output:
        tuple val(metadata), path('*.fastq.gz'), emit: reads

    script:
        if(metadata.readType == 'single') {
            """
            cat ${reads1} > ${metadata.sampleName}_R1.fastq.gz
            """
        } else if(metadata.readType == 'paired') {
            """
            cat ${reads1} > ${metadata.sampleName}_R1.fastq.gz
            cat ${reads2} > ${metadata.sampleName}_R2.fastq.gz
            """
        }
}
