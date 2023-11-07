# ngs

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

![utia-gc/ngs pipeline diagram][pipeline_diagram]

## Introduction

`ngs` is a [Nextflow](https://www.nextflow.io/) pipeline built on [utia-gc/ngs](htpps://github.com/utia-gc/ngs) for base NGS analysis.
While `ngs` can be run on any platform supported by Nextflow, it is developed for use in HPC environments and specifically [ISAAC Next Generation] at the University of Tennessee, Knoxville.

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
nextflow run utia-gc/ngs -profile exploratory
```

Once you have settled on an optimal set of parameters, rerun the pipeline without the `exploratory` profile:

```bash
nextflow run utia-gc/ngs -resume
```

Useful tip --- if a specific previous run contained the user's optimal set of parameter, or more generally if for some reason it would be advantageous to resume from some run of the pipeline other than the most recent run, then the pipeline can be resumed from any previous cached run using the RUN NAME or SESSION ID of the desired run.
Use `nextflow log` to view information about previous runs.
For example, to resume from a run named 'boring_euler':

```bash
nextflow run utia-gc/ngs -resume boring_euler
```

[pipeline_diagram]: docs/images/mermaid-diagram-2023-11-07-152502.png
