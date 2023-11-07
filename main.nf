/*
---------------------------------------------------------------------
    utia-gc/ngs
---------------------------------------------------------------------
https://github.com/utia-gc/ngs
*/

nextflow.enable.dsl=2

/*
---------------------------------------------------------------------
    RUN MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { CHECK_QUALITY  } from "./workflows/check_quality.nf"
include { GENOTYPE       } from "./workflows/genotype.nf"
include { MAP_READS      } from "./workflows/map_reads.nf"
include { PREPARE_INPUTS } from "./workflows/prepare_inputs.nf"
include { PROCESS_READS  } from "./workflows/process_reads.nf"

workflow {
    PREPARE_INPUTS(
        params.samplesheet,
        params.genome,
        params.annotations
    )
    ch_reads_raw    = PREPARE_INPUTS.out.samples
    ch_reads_raw.dump(tag: "ch_reads_raw")
    ch_genome       = PREPARE_INPUTS.out.genome
    ch_genome_index = PREPARE_INPUTS.out.genome_index
    ch_genome_dict  = PREPARE_INPUTS.out.genome_dict
    ch_annotations  = PREPARE_INPUTS.out.annotations

    PROCESS_READS(ch_reads_raw)
    ch_reads_pre_align = PROCESS_READS.out.reads_pre_align
    ch_trim_log        = PROCESS_READS.out.trim_log

    MAP_READS(
        ch_reads_pre_align,
        ch_genome,
        ch_annotations,
        params.tools.map
    )
    ch_alignmentsIndividual = MAP_READS.out.alignmentsIndividual
    ch_alignmentsMerged     = MAP_READS.out.alignmentsMerged

    // create channel with genome, genome fasta index, and genome sequence dictionary
    ch_genome
        .concat ( ch_genome_index, ch_genome_dict )
        .collect()
        .set { ch_genome_index_dict }
    GENOTYPE(
        ch_alignmentsMerged,
        ch_genome_index_dict
    )

    CHECK_QUALITY(
        ch_reads_raw,
        ch_reads_pre_align,
        ch_trim_log,
        ch_genome_index,
        ch_alignmentsIndividual,
        ch_alignmentsMerged
    )
}
