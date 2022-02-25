process AlignBowtie2 {
    tag "${metadata.sampleName}"

    container 'quay.io/biocontainers/bowtie2:2.4.5--py38hfbc8389_2'

    publishDir "${params.baseDirData}/align", mode: 'copy'

    input:
        tuple val(metadata), file(reads)
        val(runName)
        path(bt2Index)

    output:
        tuple val(metadata), path('*'), emit: sam

    script:
        if (metadata.readType == 'single') {
            argReads = "-1 ${reads}"
        } else {
            argReads = "-1 ${reads[0]} -2 ${reads[1]}"
        }

        """
        bowtie2 \
            -x ${bt2Index} \
            ${argReads} \
            -S ${runName}__bt2.sam
        """
}