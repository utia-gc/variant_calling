# Reference genomes

**IMPORTANT**: Please note that `SRAlign` **does not** package reference genome files with the workflow scripts or provide them anywhere itself. References must be acquired and supplied by the user.

## Built-in support for reference genomes

### Important information about built-in references

`SRAlign` **does not** package reference genome files with the workflow scripts or provide them anywhere itself.
*References must be acquired and supplied by the user.*
Note that this does not apply for the small test reference files which are downloaded from the remote repo [trev-f/SRAlign-test](https://github.com/trev-f/SRAlign-test) when `SRAlign` is run with `-profile test`.

### The `SRAlign` genomes config file

`SRAlign` **does** provide a [genomes config file](configs/genomes.config) for convenience.
The genomes config file consists of a single parameter, `genomes`, within the `params` scope of a Nextflow config file.
`genomes` contains key:value pairs where the key is the name of a reference genome (e.g. ce10, EB1, WBcel235) and the values are variables which set paths to the locations of specific reference files on the user's disk or a remote location.
For simplicity, currently available `genomes` key options follow the Illumina iGenomes directory and file structure conventions (see below).
This means that for many commonly used genomes, getting started is as easy as downloading the desired tar ball for that genome from the [Illumina iGenomes collection](https://support.illumina.com/sequencing/sequencing_software/igenome.html), extracting the tar ball, and supplying the correct base genome directory to `SRAlign` at run time.
See the section on default behavior below for information on how `SRAlign` uses its supplied reference genomes.

## `SRAlign` default usage

`SRAlign` accesses reference sequence and annotation files by looking for the supplied `genome` key in the `genomes` parameter and looking for a proper variable that contains the path to the desired file.

### A quick example

As an example, consider that the genomes config file contains the following (truncated) `genomes` parameter:

```Nextflow
params {
    genomes {
        'WBcel235' {
            fasta       = "${params.baseDirGenome}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
            bwa         = "${params.baseDirGenome}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BWAIndex/genome.fa"
            bowtie2     = "${params.baseDirGenome}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/genome"
        }
    }
}
```

When the user specifies WBcel235 as their desired reference genome, e.g. by supplying the command line option `--genome WBcel235`, `SRAlign` first checks that `params.genomes.WBcel235` exists.
Since this does exist in our example, `SRAlign` proceeds to run.
If `SRAlign` needs to access the fasta file within the workflow, it will use the file specified at `params.genomes.WBcel235.fasta`.
Notice that these file paths all begin with the base genome directory parameter (by default this is set to `data/references`).
If the file is found, `SRAlign` will continue to run, otherwise the command that requires the fasta file will fail and the pipeline will terminate with an error.

### Using reference genomes available in the genomes config file

This section contains step by step instructions for how to use a reference genome that is available in the [genomes config file](configs/genomes.config).

In this example, I will be aligning to the *E. coli* reference genome EB1 since it is relatively small and will serve as a quick example.

1. Create and enter a project directory for this example

    ```bash
    mkdir SRAlign_EB1_example
    cd SRAlign_EB1_example
    ```

2. Make sure the desired reference genome is available in the genomes config file. This can be done in a number of ways:
   1. See the section at the bottom of this document.
   2. Run `nextflow run trev-f/SRAlign --help` to see the available references listed on the genome parameter in the help documentation
   3. See the genome keys available in the [genomes config file](configs/genomes.config).
3. EB1 is a valid reference, so proceed. Download the reference's tar ball from the [Illumina iGenomes collection](https://support.illumina.com/sequencing/sequencing_software/igenome.html) to the directory `data/references` within your project directory.
   1. I actually prefer to have a central data directory under my home directory that keeps all of my references in a central, easy to find location, so I will download the reference files to that location instead and create a link to it in my project directory. These files could just be downloaded directly to a project `data/references` directory if that is desired.

        ```bash
        # create a data references directory within the home directory
        mkdir -p ~/data/references

        # download the EB1 tar ball to the home data references directory
        wget -P ~/data/references/ http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Escherichia_coli_K_12_DH10B_Ensembl_EB1.tar.gz

        # create a data directory within the project directory
        mkdir data

        # link project data/references to home data/references
        ln -s ~/data/references/ ./data/
        ```

   2. Extract the tar ball to the `data/references` directory

        ```bash
        tar -xvf data/references/Escherichia_coli_K_12_DH10B_Ensembl_EB1.tar.gz -C data/references/
        ```

4. Run `SRAlign` in test mode with EB1 as your genome key. This should work! Note that preseq should be skipped because the low number of aligned reads causes it to fail.

   ```bash
   nextflow run trev-f/SRAlign -profile test --genome EB1 --skipPreseq
   ```

## Illumina iGenomes

We recommend downloading common genome reference sequence and annotation files from [Illumina's iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) collection.
This collection contains reference sequences and annotations for many common model organisms downloaded and compiled from Ensembl, NCBI, or UCSC.
iGenomes contain both reference sequences and annotations in a stereotypical directory structure and file location that allows for easy programmatic access of files.
Because iGenomes are downloaded as tar balls they are easy to access and work with.
See recommendations below for adding a new genome.

`SRAlign` provides support for quasi built-in reference genomes through the provided [genomes config file](configs/genomes.config).
