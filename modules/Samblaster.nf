/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-27
Purpose: Nextflow module
*/

process Samblaster {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/samblaster:0.1.26--h9f5acd7_2'

    input:
        tuple val(metadata), stdin, val(toolIDs)

    output:
        tuple val(metadata), stdout, val(toolIDs), emit: sam

    script:
        // update toolID and set suffix
        toolIDs += 'sbl'
        suffix = toolIDs ? "__${toolIDs.join('_')}" : ''


        """
        samblaster
        """
}
