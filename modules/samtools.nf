process SAMTOOLS_SORT {
    tag "${metadata.sampleName}"

    label 'samtools'
    label 'med_mem'

    input:
        tuple val(metadata), path(bam)

    output:
        tuple val(metadata), path("*_sorted.bam"), emit: sort_star_bam

    script:
        """
        samtools sort \
            -@ ${task.cpus} \
            -o ${bam.baseName}_sorted.bam \
            ${bam} 
        """
}
