# utia_rnaseq

## Introduction

**utia_rnaseq** is a [Nextflow](https://www.nextflow.io/) pipeline for <>

## Pipeline overview

1. QC

## Quick start

### Prerequisites

1. Any POSIX compatible system (e.g. Linux, OS X, etc) with internet access

   - Run on Windows with [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/). WSL2 highly recommended.

2. [Nextflow](https://www.nextflow.io/) version >= 21.04

   - See [Nextflow Get started](https://www.nextflow.io/docs/latest/getstarted.html#) for prerequisites and instructions on installing and updating Nextflow.

3. [Singularity](https://sylabs.io)


### Get or update `utia_rnaseq`

1. Download or update `utia_rnaseq`:

    ```bash
    nextflow pull trev-f/utia_rnaseq
    ```

2. Show project info:

    ```bash
    nextflow info trev-f/utia_rnaseq
    ```

### Test `utia_rnaseq`

1. Check that `utia_rnaseq` works on your system:

      - `-profile test` uses preconfigured test parameters to run `utia_rnaseq` in full on a small test dataset stored in a remote GitHub repository.
      - Because these test files are stored in a remote repository, internet access is required to run the test.
      - For more information, see the `profiles` section of the [nextflow config file](nextflow.config) and [trev-f/utia_rnaseq-test](https://github.com/trev-f/utia_rnaseq-test).

    ```bash
    nextflow run utia_rnaseq -profile test 
    ```

### Run `utia_rnaseq`

1. Prepare the [input design csv file](docs/input_output.md).

    - Input design file must be in csv format with no whitespace.
    - Either reads (fastq or fastq.gz) or alignments (bam) are accepted.
      - If reads are supplied, can be paired or unpaired.
    - Required columns:
      - reads: lib_ID, sample_name, replicate, reads1, reads2 (optional)
      - alignments: lib_ID, sample_name, replicate, bam, tool_IDs
    - See [sample inputs](https://github.com/trev-f/utia_rnaseq-test/tree/main/inputs) in the [`utia_rnaseq-test` repository](https://github.com/trev-f/utia_rnaseq-test).
    - A template project repository can be downloaded from the [`utia_rnaseq-template` repository](https://github.com/trev-f/utia_rnaseq-template).

