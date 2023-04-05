process STAR_INDEX {
  label 'star'
  label 'sup_mem'

  //publishDir(path: "${publish_dir}/star", mode: "symlink")

  input:
  path reference
  path annotation

  output:
  path "./Star_index" , emit : star_idx

  script:
  """
  mkdir Star_index
  STAR \
  --runMode genomeGenerate \
  --genomeDir ./Star_index \
  --genomeFastaFiles ${reference} \
  --runThreadN ${task.cpus} \
  --sjdbGTFfile ${annotation} \
  --sjdbGTFtagExonParentTranscript Parent \
  --sjdbOverhang 149
  """
}

process STAR_MAP {
  label 'star'

  //publishDir(path: "${publish_dir}/star", mode: "symlink")

  input:
  tuple val(sample_id), path(reads)
  path star_idx

  output:
  tuple val(sample_id), path("*.bam") , emit : star_bam

  script:
  """
  STAR \
  --genomeDir ${star_idx} \
  --readFilesIn ${reads[0]} ${reads[1]} \
  --runThreadN ${task.cpus} \
  --outFileNamePrefix ${sample_id}_ \
  --readFilesCommand zcat \
  --limitBAMsortRAM ${task.memory.bytes} \
  --outSAMtype BAM SortedByCoordinate
  """
}


