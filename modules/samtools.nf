process BAM_TO_SAM {
  label 'samtools'
  label 'lil_mem'

  input:
  tuple val(sample_id), path(bam)

  output:
  path "*sam"

  script:
  """
  samtools view -h ${bam} > ${bam.baseName}.sam
  """
}
