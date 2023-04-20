process HTSEQ_COUNT {
  label 'htseq'
  label 'lil_mem'

  //publishDir(path: "${publish_dir}/htseq", mode: "symlink")

  input:
  tuple val(metadata), path(bam)
  path annotation

  output:
  tuple val(metadata), path("*.counts.txt"),    emit: htseq_counts

  script:
  // TODO: Find some way to incorporate strandedness
  // TODO: Add customizable options
  """
  htseq-count \
      --order=pos \
      --stranded=no \
      ${bam} \
      ${annotation} \
      > ${metadata.sampleName}.counts.txt \
      2> ${metadata.sampleName}.out
  """
}
