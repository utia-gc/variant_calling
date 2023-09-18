process bwa_mem2_mem {
    tag "${metadata.sampleName}"

    label 'bwa_mem2'

    input:
        tuple val(metadata), path(reads)
        path  index

    output:
        tuple val(metadata), file('*.sam'), emit: alignments

    script:
        // update alignments metadata
        metadata.put('format', 'SAM')
        metadata.put('sorted', false)

        // get index prefix
        def indexPrefix = index[0].toString() - ~/\.0123/

        """
        bwa-mem2 mem \
            -t ${task.cpus} \
            ${indexPrefix} \
            ${reads} \
            > ${metadata.sampleName}.sam
        """
}
