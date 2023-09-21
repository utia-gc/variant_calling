SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

ifeq ($(origin .RECIPEPREFIX), undefined)
  $(error This Make does not support .RECIPEPREFIX. Please use GNU Make 4.0 or later)
endif
.RECIPEPREFIX = >

# Run the most "up-to-date" rule
run: pipeline
.PHONY: run

pipeline: src/batch/run_ngs.bash nextflow.config data/samplesheets/samplesheet.csv
> sbatch src/batch/run_ngs.bash
.PHONY: pipeline

nextflow.config: src/nextflow/args.config src/nextflow/params.config src/nextflow/profiles.config
> touch $@

data/samplesheets/samplesheet.csv: data/reads/raw/%.fastq.gz
> bash src/bash/construct_samplesheet.bash > $@

readsdir = {{ cookiecutter.raw_reads_dir }}

data/reads/raw/%.fastq.gz :
> find \
>   $(readsdir) \
>   -maxdepth 1 \
>   -regextype egrep \
>   -regex '^/.*/.*\.fastq\.gz' \
>   -exec cp --update --preserve=all {} $(@D) \
>   \;
