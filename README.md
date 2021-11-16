# sralign
A flexible pipeline for short read alignment to a reference.

## Introduction

**sralign** is a [Nextflow](https://www.nextflow.io/) pipeline for aligning short reads to a reference. 

**sralign** is designed to be highly flexible by allowing for the easy addition of tools to the pipeline as well as serving as a starting point for genomic analyses that rely on alignment of short reads to a reference.

## Pipeline overview

1. QC of raw reads - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) & [MultiQC](https://multiqc.info/)
2. Trim raw reads - [cutadapt](https://github.com/marcelm/cutadapt)
3. Align reads - [BWA](http://bio-bwa.sourceforge.net/) -OR- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
4. Mark duplicates - [samblaster](https://github.com/GregoryFaust/samblaster)
5. QC of alignments - [Samtools](http://www.htslib.org/) & [MultiQC](https://multiqc.info/) 