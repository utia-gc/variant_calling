include { bwa_mem2_index } from '../modules/bwa_mem2_index.nf'
include { bwa_mem2_mem   } from '../modules/bwa_mem2_mem.nf'


/**
 * Workflow to run bwa-mem2.
 * 
 * Generate bwa-mem2 index and align reads to reference genome using bwa-mem2.
 * @see https://github.com/bwa-mem2/bwa-mem2
 * 
 * @take reads the reads channel of format [metadata, [R1, R2]] where R2 is optional.
 * @take genome the uncompressed reference genome sequence in fasta format.
 * @emit alignments the aligned/mapped reads channel of format [metadata, alignments] where the alignments are in unsorted SAM format and metadata has additional fields to reflect this.
 */
workflow Bwa_Mem2 {
    take:
        reads
        genome

    main:
        // build reference index
        bwa_mem2_index(genome)
        ch_index = bwa_mem2_index.out.index

        // map reads to reference
        bwa_mem2_mem(
            reads,
            ch_index
        )
        ch_alignments = bwa_mem2_mem.out.alignments

    emit:
        alignments = ch_alignments
}
