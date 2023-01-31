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
      - Allows for accessing `SRAlign` using Nextflow command by simply referring to `trev-f/SRAlign` without having to refer to the location of `SRAlign` in the system.
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
    nextflow run trev-f/SRAlign -profile test 
    ```

### Run `SRAlign`

1. Prepare the [input design csv file](docs/input_output.md).

    - Input design file must be in csv format with no whitespace.
    - Either reads (fastq or fastq.gz) or alignments (bam) are accepted.
      - If reads are supplied, can be paired or unpaired.
    - Required columns:
      - reads: lib_ID, sample_name, replicate, reads1, reads2 (optional)
      - alignments: lib_ID, sample_name, replicate, bam, tool_IDs
    - See [sample inputs](https://github.com/trev-f/SRAlign-test/tree/main/inputs) in the [`SRAlign-test` repository](https://github.com/trev-f/SRAlign-test).
    - A template project repository can be downloaded from the [`SRAlign-template` repository](https://github.com/trev-f/SRAlign-template).

2. Show all configurable options for `SRAlign` by showing a help message:

    - The most important information here is probably the list of available reference genomes.

    ```bash
    nextflow run trev-f/SRAlign --help
    ```

3. Analyze your data with `SRAlign`:

    ```bash
    nextflow run trev-f/SRAlign -profile docker --input <input.csv> --genome <valid genome key>
    ```

## Tips for running Nextflow and `SRAlign`

`SRAlign` is designed to be highly configurable, meaning that its default behavior can be changed by supplying any of a number of configurable parameters. These can be [supplied in a number of ways](https://www.nextflow.io/docs/latest/config.html#configuration-file) that have a specific hierarchy of precedence.

- Show configurable parameters by showing command line help documentation: `nextflow run trev-f/SRAlign --help`
- Nextflow arguments always begin with a single dash, e.g. `-profile`.
- Pipeline parameters specified at the command line always begin with a double dash, e.g. `--input`.
  - Parameters specified at the command line always have the highest precedence. They will overwrite parameters specified in any config or params files.
  - I recommend specifying required parameters (i.e. `--input` and `--genome`) and up to a few others at the command line in this manner. Specifying more than this at the command line gets unwieldy.
- A custom config or parameters file is a good option for cases where you want to supply more parameters than can comfortably be done at the command line or you want to use the same custom parameters in multiple runs.
  - For a config file, use the [params scope](https://www.nextflow.io/docs/latest/config.html#scope-params)
  - For a JSON/YAML parameters file, see the [Nextflow CLI docs](https://www.nextflow.io/docs/latest/cli.html?highlight=params%20file#run).

## Additional documentation

Additional documentation can be found in [docs](docs/).

Quick links:

- [Information about SRAlign](docs/SRAlign_info.md)
- [Detailed usage](docs/usage.md)
- [Inputs and outputs](docs/input_output.md)
- [Pipeline configuration](docs/configuration.md)
- [Reference genomes](docs/reference_genomes.md)
