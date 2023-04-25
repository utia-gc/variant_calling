process STAR_INDEX {
    label 'star'
    label 'big_mem'
  
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
        --sjdbOverhang ${params.sjdbOverhang}
        """
}

process STAR_MAP {
      label 'star'
      label 'sup_mem'
    
      container = "quay.io/biocontainers/star:2.7.10b--h9ee0642_0"
    
      //publishDir(path: "${publish_dir}/star", mode: "symlink")
    
      input:
          tuple val(metadata), path(reads)
          path star_idx
    
      output:
          tuple val(metadata), path("*.bam") , emit : star_bam
    
      script:
          """
          STAR \
          --genomeDir ${star_idx} \
          --readFilesIn ${reads} \
          --runThreadN ${task.cpus} \
          --outFileNamePrefix ${metadata.sampleName}_ \
          --readFilesCommand zcat \
          --outSAMtype BAM Unsorted
          """
}


