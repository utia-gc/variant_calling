process Fastp {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0'

    input:
        tuple val(metadata), file(reads)

    output:
        tuple val(metadata), path('*.fastq.gz'), emit: trimReads

    script:
        if (metadata.readType == 'single') {
            adapterTrim = params.adapterR1 ? "--adapter_sequence ${params.adapterR1}" : ''

            """
            fastp \
                ${adapterTrim} \
                -i ${reads} \
                -o ${metadata.sampleName}__fsp_R1.fastq.gz
            """
        } else {
            adapterTrimR1 = params.adapterR1 ? "--adapter_sequence ${params.adapterR1}" : ''
            adapterTrimR2 = params.adapterR2 ? "--adapter_sequence_r2 ${params.adapterR2}" : ''

            """
            fastp \
                ${adapterTrimR1} \
                ${adapterTrimR2} \
                -i ${reads[0]} \
                -I ${reads[1]} \
                -o ${metadata.sampleName}__fsp_R1.fastq.gz \
                -O ${metadata.sampleName}__fsp_R2.fastq.gz
            """
        }
}