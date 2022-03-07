/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-03-06
Purpose: BAlign reads to genome with hisat2
*/

process Hisat2Align {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/hisat2:2.2.1--h87f3376_4'

    label 'cpu_mid'

    input:
        tuple val(metadata), file(reads), val(toolIDs)
        path ht2Indexes

    output:
        tuple val(metadata), stdout, val(toolIDs), emit: sam

    script:
        // set reads arguments
        if (metadata.readType == 'single') {
            argReads = "-U ${reads}"
        } else {
            argReads = "-1 ${reads[0]} -2 ${reads[1]}"
        }

        // update toolID and set suffix
        toolIDs += "ht2-${params.genome}"
        suffix = toolIDs ? "__${toolIDs.join('_')}" : ''

        // set index base name
        def ht2IndexBaseName = ht2Indexes[0].toString() - ~/.rev.\d.ht2?/ - ~/.\d.ht2?/ - ~/.fa?/

        // set arguments
        def options = task.ext.args ?: ''

        """
        hisat2 \
            ${options} \
            -x ${ht2IndexBaseName} \
            ${argReads}
        """
}
