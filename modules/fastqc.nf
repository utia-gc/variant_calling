process FASTQC {
  label 'fastqc'
  label 'lil_mem'

  container = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

  publishDir(path: "${publish_dir}/raw_qc", mode: "symlink")

  input:
  tuple val(sample_id), path(reads)

  output:
  path("fastqc_${sample_id}_logs"), emit: fastq_ch

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -q ${reads}
  """
}
