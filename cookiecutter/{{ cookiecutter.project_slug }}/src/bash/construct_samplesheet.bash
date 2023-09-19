#!/usr/bin/env bash

# write header for paired-end reads
echo "sampleName,reads1,reads2"

# add sample name and path to fastq files
for FQ in $(ls data/reads/raw/*.fastq.gz | sed 's/_R[12]_001.fastq.gz//' | uniq);
do
    BASE=$(basename $FQ)
    echo "${BASE},${FQ}_R1_001.fastq.gz,${FQ}_R2_001.fastq.gz"
done
