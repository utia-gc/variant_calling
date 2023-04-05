process MULTIQC {
  label 'multiqc'
  label 'lil_mem'

  publishDir(path: "${publish_dir}/raw_qc", mode: "symlink")

  input:
  path('*')

  script:
  """
  multiqc .
  """
}
