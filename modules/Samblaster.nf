/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-27
Purpose: Nextflow module
*/

process Samblaster {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/samblaster:0.1.26--h9f5acd7_2'

    input:
        tuple val(metadata), stdin

    output:
        tuple val(metadata), stdout, emit: sam

    script:
        """
        samblaster
        """
}
