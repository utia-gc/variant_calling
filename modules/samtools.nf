process SAMTOOLS_SORT {
    label 'samtools'
    label 'med_mem'
 
    container = "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
 
    input:
        tuple val(metadata), path(bam)
  
    output:
        tuple val(metadata), path("*sorted.bam") , emit : sort_star_bam
  
    script:
        """
        samtools sort -@ ${task.cpus} ${bam} > ${bam.baseName}.sorted.bam
        """
}
