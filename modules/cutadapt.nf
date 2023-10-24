process cutadapt {
    tag "${metadata.sampleName}"
    
    label 'cutadapt'
    label 'lil_mem'

    publishDir(
        path:    "${params.publishDirReports}/reads/trim",
        mode:    "${params.publishMode}",
        pattern: '*_cutadapt-log.txt'
    )

    input:
        tuple val(metadata), path(reads1), path(reads2)
        val r1_adapter
        val r2_adapter
        val minimum_length

    output:
        tuple val(metadata), path("*_trimmed_R1.fastq.gz"), path("*_trimmed_R2{.fastq.gz,.NOFILE}"), emit: reads
        path("*_cutadapt-log.txt"), emit: log

    script:
        String stemName = MetadataUtils.buildStemName(metadata)
        String reads1NewName = "${stemName}_${metadata.trimStatus}_R1.fastq.gz"

        if(metadata.readType == 'single') {
            """
            mv ${reads1} ${reads1NewName}

            cutadapt \
                -a ${r1_adapter} \
                -m ${minimum_length} \
                -o ${stemName}_trimmed_R1.fastq.gz \
                ${reads1NewName} \
                > ${stemName}_cutadapt-log.txt

            cp ${reads2} ${stemName}_trimmed_R2.NOFILE
            """
        } else if(metadata.readType == 'paired') {
            String reads2NewName = "${stemName}_${metadata.trimStatus}_R2.fastq.gz"

            """
            mv ${reads1} ${reads1NewName}
            mv ${reads2} ${reads2NewName}

            cutadapt \
                -a ${r1_adapter} \
                -A ${r2_adapter} \
                -m ${minimum_length} \
                -o ${stemName}_trimmed_R1.fastq.gz \
                -p ${stemName}_trimmed_R2.fastq.gz \
                ${reads1NewName} ${reads2NewName} \
                > ${stemName}_cutadapt-log.txt
            """
        }
}
