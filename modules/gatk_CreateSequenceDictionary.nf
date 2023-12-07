/**
 * Create reference sequence dictionary file
 *
 * Create the ref.dict file from reference fasta sequence that is required for many GCTA tools.
 * @see https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
 * 
 * @input genome the uncompressed reference genome sequence in fasta format.
 * @emit dict the reference sequence .dict file
 */
process gatk_CreateSequenceDictionary {
    tag "${genome.name}"

    label 'gatk'

    label 'med_mem'

    input:
        path genome

    output:
        path '*.dict', emit: dict

    script:
        String args = new Args(task.ext).buildArgsString()

        """
        gatk CreateSequenceDictionary \\
            --REFERENCE ${genome} \\
            ${args}
        """
}
