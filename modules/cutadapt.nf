process CUTADAPT {
    tag "${metadata.sampleName}"
    
    label 'cutadapt'
    label 'lil_mem'

    publishDir(path: "${publish_dir}/cutadapt", mode: "symlink")

    input:
        tuple val(metadata), path(reads)
        val r1_adapter
        val r2_adapter
        val minimum_length

    output:
        tuple val(metadata), path("${metadata.sampleName}_*.fastq.gz"), emit: reads
        path("${metadata.sampleName}_cutadapt-log.txt"), emit: log

    script:
        if (metadata.readType == 'single') {
            """
            cutadapt \
                -a ${r1_adapter} \
                -m ${minimum_length} \
                -o ${metadata.sampleName}_R1.fastq.gz \
                ${reads} \
                > ${metadata.sampleName}_cutadapt-log.txt
            """
        } else if(metadata.readType == 'paired') {
            """
            cutadapt \
                -a ${r1_adapter} \
                -A ${r2_adapter} \
                -m ${minimum_length} \
                -o ${metadata.sampleName}_R1.fastq.gz \
                -p ${metadata.sampleName}_R2.fastq.gz \
                ${reads} \
                > ${metadata.sampleName}_cutadapt-log.txt
            """
        }
}
