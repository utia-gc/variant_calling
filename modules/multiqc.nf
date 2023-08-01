process MULTIQC {
  label 'multiqc'
  label 'lil_mem'

  publishDir(path: "${publish_dir}/qc/${fileName}", mode: "copy")

  input:
      path('*')
      val fileName

  output:
      path("*html"),          emit: report
      path("multiqc_data/*"), emit: data

  script:
      """
      multiqc \
        --filename ${fileName} \
        .
      """
}
