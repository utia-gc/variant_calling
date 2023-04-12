process MULTIQC {
  label 'multiqc'
  label 'lil_mem'

  container = "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"  

  publishDir(path: "${publish_dir}/qc/${outdir_name}", mode: "copy")

  input:
      path('*')
      val outdir_name

  output:
      path("*html")

  script:
      """
      multiqc .
      """
}
