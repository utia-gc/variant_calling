process AlignBowtie2 {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/bowtie2:2.4.5--py38hfbc8389_2'

    publishDir "${params.baseDirData}/align", mode: 'copy'

    input:
        tuple val(metadata), file(reads)
        val runName
        path bt2Indexes

    output:
        tuple val(metadata), stdout, emit: sam

    script:
        // set reads arguments
        if (metadata.readType == 'single') {
            argReads = "-1 ${reads}"
        } else {
            argReads = "-1 ${reads[0]} -2 ${reads[1]}"
        }


        // determine index base name
        bt2IndexBaseName = bt2Indexes[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/ - ~/.fa?/

        """
        bowtie2 \
            -x ${bt2IndexBaseName} \
            ${argReads}
        """
}
