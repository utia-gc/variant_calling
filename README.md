# NGS

## Introduction

`ngs` is a [Nextflow](https://www.nextflow.io/) pipeline for base NGS analysis.
`ngs` is primarily intended to be used as a starting point for building more specific NGS analysis pipelines, e.g. for RNA-seq, WGS genotyping, snATAC-seq, etc.
While `ngs` can be run on any platform supported by Nextflow, it is developed for use in HPC environments and specifically [ISAAC Next Generation] at the University of Tennessee, Knoxville.

## Pipeline overview

1. Trim/filter reads
2. Reads QC
3. Map reads
4. Process alignments
5. Full workflow QC

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
   nextflow run utia-gc/ngs -profile nf_test 
   ```

### Run `ngs`

TODO
