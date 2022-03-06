/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-03-06
Purpose: Build bowtie2 index
*/

process Bowtie2Build {
    tag "${params.genome}"

    container 'quay.io/biocontainers/bowtie2:2.4.5--py38hfbc8389_2'

    input:
        path reference

    output:
        path '*', emit: bowtie2Index

    script:
        """
        bowtie2-build ${reference}
        """
}
