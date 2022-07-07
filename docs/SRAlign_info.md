# SRAlign Info

## Why SRAlign?

`SRAlign` is designed to be a highly flexible, highly configurable, and reproducible general purpose next-generation sequencing pipeline.

I designed `SRAlign` when I got tired of manually running the same commands over and over again for a variety of NGS datasets: FastQC to check sequencing reads, adapter trimming, FastQC again on those trimmed reads, aligning to a reference genome, marking duplicate alignments, sorting and indexing bam files, checking alignment stats.
This process is incredibly inefficient on a number of levels, including how mistake prone the process is to mistakes and the amount of time that has to be devoted to entering commands and waiting for them to run.
Further, this process is essentially the same regardless of the type of datasets I was analyzing from ChIP-seq and ATAC-seq to RNA-seq.

I have two main motivations in mind with `SRAlign`. 1) To make a pipeline that performs the central tasks of aligning short reads to a reference genome along with the necessary pre- and postprocessing steps and with extensive QC reporting. 2) To make a pipeline that can be easily edited to add new functionalities and, where appropriate, extended into a new pipeline.
