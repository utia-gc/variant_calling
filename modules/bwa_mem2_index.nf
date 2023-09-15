process bwa_mem2_index {
    tag "${in_fasta.name}"

    label 'bwa_mem2'

    input:
        path in_fasta

    output:
        path "${in_fasta.name}*", emit: index

    script:
        """
        bwa-mem2 index ${in_fasta}
        """
}
