# ngs

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

`ngs` is a [Nextflow](https://www.nextflow.io/) pipeline built on [utia-gc/ngs](htpps://github.com/utia-gc/ngs) for base NGS analysis.
While `ngs` can be run on any platform supported by Nextflow, it is developed for use in HPC environments and specifically [ISAAC Next Generation] at the University of Tennessee, Knoxville.

> [!WARNING]
> Do not directly run `utia-gc/ngs`. This pipeline is designed to be a starting point for pipelines dedicated to specific analyses. It is not meant to be run itself. It does not have any scheduled or tagged releases, and as such it cannot be used reproducibly.

### Pipeline overview

```mermaid
flowchart LR
    %% list all the input files
    samplesheet>"Samplesheet"]
    genome_fasta>"
        Genome
        FASTA
    "]
    annotations_gtf>"
        Annotations
        GTF
    "]

    %% list all the internal Nextflow channels
    raw_reads[("
        Raw
        reads
    ")]
    prealign_reads[("
        Prealign
        reads
    ")]
    trim_log[("
        Trim
        log
    ")]
    individual_alignments[("
        Individual
        alignments
    ")]
    merged_alignments[("
        Merged
        alignments
    ")]

    %% list all the Nextflow processes
    fastp{"fastp"}
    cutadapt{"cutadapt"}
    fastqc{"FastQC"}
    seq_depth{"
        Sequencing
        Depth
    "}
    bwa_mem2{"bwa-mem2"}
    STAR{"STAR"}
    samtools_sort{"
        samtools
        sort
        index
    "}
    gatk_MergeSamFiles{"
        gatk
        MergeSamFiles
    "}
    gatk_MarkDuplicates{"
        gatk
        MarkDuplicates
    "}
    samtools_idxstats{"
        samtools
        idxstats
    "}
    samtools_flagstat{"
        samtools
        flagstat
    "}
    samtools_stats{"
        samtools
        stats
    "}

    %% list all subgraphs for Nextflow subworkflows/workflows with options
    subgraph inputs["Input Files"]
    samplesheet
    genome_fasta
    annotations_gtf
    end
    subgraph trim_reads["Trim Reads"]
    fastp
    cutadapt
    end
    subgraph map_reads["Map Reads"]
    bwa_mem2
    STAR
    end
    subgraph publish_reports["Publish Reports"]
    reads_mqc
    alignments_mqc
    full_mqc
    end
    subgraph publish_data["Publish Data"]
    alignments
    end

    %% list all the published reports files
    reads_mqc((("
        Reads
        MultiQC
    ")))
    alignments_mqc((("
        Alignments
        MultiQC
    ")))
    full_mqc((("
        Full MultiQC
    ")))

    %% list all the published data files
    alignments[["
        Alignments
    "]]

    %% reads processing workflow
    samplesheet --> raw_reads
    raw_reads --- trim_reads --> prealign_reads

    %% reads QC workflow
    raw_reads --- fastqc --x reads_mqc
    prealign_reads --- fastqc --x reads_mqc
    trim_reads --> trim_log --x reads_mqc
    raw_reads --- seq_depth --x reads_mqc
    prealign_reads --- seq_depth --x reads_mqc

    %% reads mapping workflow
    genome_fasta --- map_reads
    annotations_gtf --- map_reads
    prealign_reads --- map_reads

    %% alignments processing workflow
    map_reads --- samtools_sort --> individual_alignments
    individual_alignments --- gatk_MergeSamFiles --- gatk_MarkDuplicates --> merged_alignments
    merged_alignments --x alignments

    %% alignments QC workflow
    individual_alignments --- samtools_idxstats --x alignments_mqc
    individual_alignments --- samtools_flagstat --x alignments_mqc
    merged_alignments --- samtools_stats --x alignments_mqc

    %% Full MultiQC
    reads_mqc --x full_mqc
    alignments_mqc --x full_mqc
```

## Quick start

### Prerequisites

1. Any POSIX compatible system (e.g. Linux, OS X, etc) with internet access

   - Run on Windows with [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/). WSL2 highly recommended.

2. [Nextflow](https://www.nextflow.io/) version >= 21.04

   - See [Nextflow Get started](https://www.nextflow.io/docs/latest/getstarted.html#) for prerequisites and instructions on installing and updating Nextflow.

3. [Singularity](https://sylabs.io)

### Get or update `ngs`

1. Download or update `ngs`:

    ```bash
    nextflow pull utia-gc/ngs
    ```

2. Show project info:

    ```bash
    nextflow info utia-gc/ngs
    ```

### Test `ngs`

1. Check that `ngs` works on your system:

   - `-profile nf_test` uses preconfigured test parameters to run `ngs` in full on a small test dataset stored in a remote GitHub repository.
   - Because these test files are stored in a remote repository, internet access is required to run the test.
   - For more information, see the `profiles` section of the [nextflow config file](nextflow.config).

   ```bash
   nextflow run utia-gc/ngs \
      -revision main \
      -profile nf_test
   ```

   > [!IMPORTANT]
   > In accordance with best practices for reproducible analysis, always use the `-revision` option in `nextflow run` to specify a tagged and/or released version of the pipeline.

### Run `ngs`

TODO

## Recommended usage

### The `exploratory` profile

Analysis of single cell RNA-seq data frequently requires an exploratory stage that involves iterating through various parameter options and scrutinizing their effects on important QC metrics before settling on a final set of parameters.
To help facilitate this crucial process, we have included an `exploratory` profile option which implements the following features:

- Data and reports are published within time-stamped subdirectories of `exploratory` structured as follows: `<current working directory>/exploratory/<timestamp>_<project title>`. This allows the user to see a chronological log of their changes and gives the option to put a brief description of changes in the project title.
- Results are published as symbolic links as opposed to the default behavior of copying published results. This prevents the user's working directory from being bloated with duplicates of data.
- Sets the `-resume` flag in Nextflow through the profile so that it does not need to be supplied at the command line. This allows for faster iteration and exploration as results from intensive processes are used from their cached location instead of being reproduced.
For more info on Nextflow's resume feature, checkout these articles on [demistifying](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html) and [troubleshooting](https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html) Nextflow resume.

Once the user finishes exploring and has decided on a final set of parameters, those parameters should be specified during an explicitly resumed run of the pipeline without the `exploratory` profile.
By default this will rerun the pipeline and publish results by copying them into the user's specified data and report publishing directories (see [output documentation](docs/output.md)).
This serves the dual purpose of saving time by not repeating logged tasks while aiding in data persistence.

#### Example usage

During exploratory analysis, iteratively make changes to parameters and run the pipeline with the `exploratory` profile:

```bash
nextflow run utia-gc/ngs \
   -revision main \
   -profile exploratory
```

Once you have settled on an optimal set of parameters, rerun the pipeline without the `exploratory` profile:

```bash
nextflow run utia-gc/ngs \
   -revision main \
   -resume
```

Useful tip --- if a specific previous run contained the user's optimal set of parameter, or more generally if for some reason it would be advantageous to resume from some run of the pipeline other than the most recent run, then the pipeline can be resumed from any previous cached run using the RUN NAME or SESSION ID of the desired run.
Use `nextflow log` to view information about previous runs.
For example, to resume from a run named 'boring_euler':

```bash
nextflow run utia-gc/ngs \
   -revision main \
   -resume boring_euler
```
