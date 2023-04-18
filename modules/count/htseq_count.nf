process HTSEQ_COUNT {
  label 'htseq'
  label 'lil_mem'

  //publishDir(path: "${publish_dir}/htseq", mode: "symlink")

  input:
  tuple val(sample_id), path(bam)
  path annotation

  output:
  path "*.counts.txt" , emit : htseq_counts

  script:
  """
  htseq-count \
      --format=bam \
      --order=pos \
      --stranded=no \
      --type=gene \
      --idattr=ID \
      ${bam} \
      ${annotation} \
      > ${sample_id}.counts.txt \
      2> ${sample_id}.out
  """
}
