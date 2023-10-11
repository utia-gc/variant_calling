/**
 * Process to compute average sequencing depth.
 *
 * Compute average sequencing depth of fastq files for a reference genome.
 * 
 * @input reads the reads channel of format [metadata, [R1, R2]] where R2 is optional.
 * @input fai the fasta index of the reference genome.
 * @emit depth path to file of average sequencing depth.
 */
process sequencing_depth {
    tag "${metadata.sampleName}"

    label 'base'

    input:
        tuple val(metadata), path(reads1), path(reads2)
        path fai

    output:
        path '*_seq-depth.txt', emit: depth

    script:
        String stemName = metadata.lane ? "${metadata.sampleName}_S${metadata.sampleNumber}_L${metadata.lane}" : "${metadata.sampleName}"
        def reads = (metadata.readType == 'single') ? "${reads1}" : "${reads1} ${reads2}"

        """
        # initialize array of reads
        READS=(${reads})

        # count number bases in fastqs
        n_reads_bases=0
        for read in \$READS; do
            n_read_bases=\$(zcat \$read | awk 'NR % 4 == 2{LEN+=length(\$0)}END{print LEN}')
            n_reads_bases=\$((n_reads_bases + n_read_bases))
        done

        # compute length of genome
        n_bases_haploid_genome=\$(awk '{s+=\$2}END{print s}' ${fai})

        # compute sequencing depth
        depth=\$(awk \
            -v n_reads_bases=\$n_reads_bases \
            -v n_genome_bases=\$n_bases_haploid_genome \
            'BEGIN { print( n_reads_bases / n_genome_bases )}')

        # write file
        echo -e "Sequencing_depth\t\$depth" > ${stemName}_seq-depth.txt
        """
}
