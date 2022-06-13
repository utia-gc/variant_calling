# SRAlign
A flexible pipeline for short read alignment to a reference with extensive QC reporting.

## Introduction

**SRAlign** is a [Nextflow](https://www.nextflow.io/) pipeline for aligning short reads to a reference. 

**SRAlign** is designed to be highly flexible by allowing for the easy addition of tools to the pipeline as well as serving as a starting point for genomic analyses that rely on alignment of short reads to a reference.

## Pipeline overview

1. QC of raw reads - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) & [MultiQC](https://multiqc.info/)
2. Trim raw reads - [cutadapt](https://github.com/marcelm/cutadapt)
3. Align reads - [BWA](http://bio-bwa.sourceforge.net/) -OR- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
4. Mark duplicates - [samblaster](https://github.com/GregoryFaust/samblaster)
5. QC of alignments - [Samtools](http://www.htslib.org/) & [MultiQC](https://multiqc.info/) 

## Quick start

1. [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
2. [Install Docker](https://docs.docker.com/engine/install/)
3. Download **SRAlign**:
    ```
    git clone https://github.com/trev-f/SRAlign.git
    ```
4. Run **SRAlign** in test mode:
    ```
    nextflow run SRAlign -profile test 
    ```
5. Run your analysis:
    ```
    nextflow run SRAlign -profile docker --input <input.csv> --genome <valid genome key>
    ```

Detailed documentation can be found in [docs](docs/) and [usage](docs/usage.md)
