include { star_genomeGenerate } from '../modules/star_genomeGenerate.nf'
include { star_runMapping     } from '../modules/star_runMapping.nf'

/**
 * Workflow to run STAR.
 * 
 * Generate STAR index and align reads to reference genome using STAR.
 * @see https://github.com/alexdobin/STAR
 * 
 * @take reads the reads channel of format [metadata, [R1, R2]] where R2 is optional.
 * @take genome the uncompressed reference genome sequence in fasta format.
 * @take annotationsGTF the uncompressed reference annotations in GTF format.
 * @emit alignments the aligned/mapped reads channel of format [metadata, alignments] where the alignments are in unsorted SAM format.
 */

workflow Star {
    take:
        reads
        genome
        annotationsGTF

    main:
        // build reference index
        star_genomeGenerate(
            genome,
            annotationsGTF
        )

        // map reads to reference
        star_runMapping(
            reads,
            star_genomeGenerate.out.index
        )

    emit:
        alignments = star_runMapping.out.alignments
}
