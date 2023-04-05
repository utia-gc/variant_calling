process DESEQ2_FROM_HTSEQ {
  label 'deseq2'
  label 'lil_mem'

  publishDir(path: "${publish_dir}/deseq2", mode: "symlink")

  input:
  path('*')
  path metadata

  script:
  """
  #!/usr/bin/env Rscript
  library(DESeq2, quietly = T)

  df = read.csv("${metadata}", header=T)
  print(df)
  dds = DESeqDataSetFromHTSeqCount(sampleTable=df, design = ~ condition)
  rld <- rlog(dds, blind = FALSE)

  pdf(file="pca_plot.pdf")
  print(plotPCA(rld, intgroup = c("condition")))
  dev.off()
  
  dds = DESeq(dds)
  dds_res = results(dds, alpha = 0.05, contrast = c("condition", "max2", "control"))

  pdf(file="ma_plot.pdf")
  print(plotMA(dds_res, ylim=c(-6,6)))
  dev.off()

  dds_res_sig = as.data.frame(dds_res[which(dds_res\$padj < 0.05),])
  dds_res_sig = dds_res_sig[order(dds_res_sig\$padj),]
  """
}
