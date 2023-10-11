process fastp {
    tag "${metadata.sampleName}"

    label 'fastp'

    label 'med_cpu'
    label 'med_mem'

    publishDir(
        path:    "${params.publishDirReports}/reads/trim",
        mode:    "${params.publishMode}",
        pattern: '*_fastp-log.json'
    )

    input:
        tuple val(metadata), path(reads1), path(reads2)

    output:
        tuple val(metadata), path("*_trimmed_R1.fastq.gz"), path("*_trimmed_R2{.fastq.gz,.NOFILE}"), emit: reads
        path("*_fastp-log.json"), emit: log

    script:
        if(metadata.readType == 'single') {
            """
            fastp \
                --thread ${task.cpus} \
                --in1 ${reads1} \
                --out1 ${metadata.sampleName}_trimmed_R1.fastq.gz \
                --json ${metadata.sampleName}_fastp-log.json

            cp ${reads2} ${metadata.sampleName}_trimmed_R2.NOFILE
            """
        } else if(metadata.readType == 'paired') {
            """
            fastp \
                --thread ${task.cpus} \
                --in1 ${reads1} \
                --in2 ${reads2} \
                --out1 ${metadata.sampleName}_trimmed_R1.fastq.gz \
                --out2 ${metadata.sampleName}_trimmed_R2.fastq.gz \
                --json ${metadata.sampleName}_fastp-log.json
            """
        }
}
