process samtools_faidx {
    tag "${fasta}"

    label 'samtools'

    input:
        path fasta

    output:
        path '*.fai', emit: fai

    script:
        """
        samtools faidx ${fasta}
        """
}
