process fastp {
    tag "${metadata.sampleName}"

    label 'fastp'
    label 'lil_mem'

    input:
        tuple val(metadata), path(reads)

    output:
        tuple val(metadata), path("*_trimmed_R{1,2}.fastq.gz"), emit: reads
        path("*_fastp-log.json"), emit: log

    script:
        if(metadata.readType == 'single') {
            """
            fastp \
                --thread ${task.cpus} \
                --in1 ${reads[0]} \
                --out1 ${metadata.sampleName}_trimmed_R1.fastq.gz \
                --json ${metadata.sampleName}_fastp-log.json
            """
        }
}
