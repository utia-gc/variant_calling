/**
 * Compress, sort, and/or index alignments (SAM or BAM) with samtools.
 *
 * @input tuple the aligned/mapped reads channel of format [metadata, alignments] where the alignments are in sorted or unsorted SAM or BAM format and metadata has fields to reflect this.
 * @output bamIndexed the sorted and indexed aligned/mapped reads channel of format [metadata, BAM, BAM.BAI] where the alignments are in sorted BAM format and metadata has fields to reflect this.
 */
process samtools_sort_index {
    tag "${metadata.sampleName}"

    label 'samtools'

    publishDir(
        path:    "${params.publishDirData}/alignments",
        mode:    "${params.publishMode}",
        pattern: "${metadata.sampleName}_sorted.bam*"
    )

    input:
        tuple val(metadata), path(alignments)

    output:
        tuple val(metadata), path("${metadata.sampleName}_sorted.bam"), path("${metadata.sampleName}_sorted.bam.bai"), emit: bamIndexed

    script:
        if(!metadata.sorted) {
            // update alignments metadata
            metadata.put('format', 'BAM')
            metadata.put('sorted', true)
            // sort alignments if not sorted, then index BAM.
            // writes a BAM regardless of input format.
            """
            samtools sort \
                -@ ${task.cpus} \
                -O bam \
                -o ${metadata.sampleName}_sorted.bam \
                ${alignments}

            samtools index \
                -@ ${task.cpus} \
                ${metadata.sampleName}_sorted.bam \
            """
        } else if(metadata.format == 'BAM' && metadata.sorted) {
            // index BAM
            """
            # produce file with stereotypical name for output
            mv -f ${alignments} ${metadata.sampleName}_sorted.bam

            samtools index \
                -@ ${task.cpus} \
                ${metadata.sampleName}_sorted.bam
            """
        }
}
