# SRAlign

A flexible pipeline for short read alignment to a reference with extensive QC reporting.

## Introduction

**SRAlign** is a [Nextflow](https://www.nextflow.io/) pipeline for aligning short reads to a reference.

**SRAlign** is designed to be highly flexible by allowing for the easy addition of tools to the pipeline as well as serving as a starting point for genomic analyses that rely on alignment of short reads to a reference.

## Pipeline overview

1. Trim reads
2. QC of reads
   1. Raw reads FastQC
   2. Trim reads FastQC
   3. Summary MultiQC
3. Align reads
    1. Align to reference genome/transcriptome
    2. Check contamination
4. Preprocess alignments
   1. Mark duplicates
   2. Compress sam to bam
   3. Index bam
5. QC of alignments
   1. samtools stats
   2. Samtools index stats
   3. Percent duplicates
   4. Percent aligned to contamination reference
   5. Summary MultiQC
6. Library complexity and reproducibility
   1. Preseq library complexity
   2. DeepTools correlation
   3. DeepTools PCA
7. Full pipeline MultiQC

## Quick start

### Prerequisites

1. Any POSIX compatible system (e.g. Linux, OS X, etc) with internet access

   - Run on Windows with [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/). WSL2 highly recommended.

2. [Nextflow](https://www.nextflow.io/) version >= 21.04

   - See [Nextflow Get started](https://www.nextflow.io/docs/latest/getstarted.html#) for prerequisites and instructions on installing and updating Nextflow.

3. [Docker](https://docs.docker.com/)

    - I recommend Docker Desktop for OS X or Windows users

### Get or update `SRAlign`

1. Download or update `SRAlign`:

    - Downloads the project into `$HOME/.nextflow/assets`
    - Useful for quickly downloading and easily running a project.
      - To customize or expand `SRAlign`, see the documentation on [customizing or expanding `SRAlign`](docs/customize_expand.md).

    ```bash
    nextflow pull trev-f/SRAlign
    ```

2. Show project info:

    ```bash
    nextflow info trev-f/SRAlign
    ```

### Test `SRAlign`

1. Check that `SRAlign` works on your system:

    - `-profile test` uses preconfigured test parameters to run `SRAlign` in full on a small test dataset stored in a remote GitHub repository.
      - Because these test files are stored in a remote repository, internet access is required to run the test.
      - For more information, see the `profiles` section of the [nextflow config file](nextflow.config) and [trev-f/SRAlign-test](https://github.com/trev-f/SRAlign-test).

    ```bash
    nextflow run SRAlign -profile test 
    ```


2. Run your analysis:
    ```
    nextflow run SRAlign -profile docker --input <input.csv> --genome <valid genome key>
    ```

Detailed documentation can be found in [docs](docs/) and [usage](docs/usage.md)
