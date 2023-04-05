process CUTADAPT_ADAPTERS {
  label 'cutadapt'
  label 'lil_mem'

  //publishDir(path: "${publish_dir}/cutadapt", mode: "symlink")

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("cut*") , emit : reads

  script:
  forward = "cut_${reads[0]}"
  reverse = "cut_${reads[1]}"
  """
  cutadapt \
    -a ${params.r1_adapter} \
    -A ${params.r2_adapter} \
    -m 30 \
    -o $forward \
    -p $reverse \
    $reads 
  """
}
