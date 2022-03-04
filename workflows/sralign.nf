/*
=====================================================================
    HANDLE INPUTS
=====================================================================
*/

/*
    ---------------------------------------------------------------------
    Design and Inputs
    ---------------------------------------------------------------------
*/

// check design file
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input design file not specified!'
}


// set input design name
inName = params.input.take(params.input.lastIndexOf('.')).split('/')[-1]


/*
    ---------------------------------------------------------------------
    Genomes and References
    ---------------------------------------------------------------------
*/

// check genome
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "Genome '${params.genome}' is not available.\n\tAvailable genomes include: ${params.genomes.keySet().join(", ")}"
}


// set genome variables
bt2Index = params.genome ? params.genomes[params.genome].bowtie2 ?: false : false


/*
=====================================================================
    MAIN WORKFLOW
=====================================================================
*/

include { ParseDesignSWF  as ParseDesign  } from '../subworkflows/ParseDesignSWF.nf'
include { RawReadsQCSWF   as RawReadsQC   } from '../subworkflows/RawReadsQCSWF.nf'
include { TrimReadsSWF    as TrimReads    } from '../subworkflows/TrimReadsSWF.nf'
include { TrimReadsQCSWF  as TrimReadsQC  } from '../subworkflows/TrimReadsQCSWF.nf'
include { AlignBowtie2SWF as AlignBowtie2 } from '../subworkflows/AlignBowtie2SWF.nf'
include { SamStatsQCSWF   as SamStatsQC   } from '../subworkflows/SamStatsQCSWF.nf'


workflow sralign {

    /*
    ---------------------------------------------------------------------
        Read design file, parse sample names and identifiers, and stage reads files
    ---------------------------------------------------------------------
    */

    // run subworkflow: Parse design file
    ParseDesign(ch_input)

    // collect output: Raw reads
    ch_rawReads = ParseDesign.out.rawReads


    /*
    ---------------------------------------------------------------------
        Raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipRawFastQC) {
        // Subworkflow: Raw reads fastqc and mulitqc
        RawReadsQC(ch_rawReads, inName)
    }


    /*
    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipTrimReads) {
        // Subworkflow: Trim raw reads
        TrimReads(ch_rawReads)
        ch_trimReads = TrimReads.out.trimReads

        if (!params.skipTrimReadsQC) {
            // Subworkflow: Trimmed reads fastqc and multiqc
            TrimReadsQC(ch_trimReads, inName) 
        }
    } 


    /*
    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */

    // run subworkflow: Align trimmed reads to genome, mark dups, sort and compress sam, and index bam
    AlignBowtie2(ch_trimReads, bt2Index)

    // collect output: Indexed bams
    ch_indexedBam = AlignBowtie2.out.bamBai

    // run subworkflow: Samtools stats and samtools idxstats and multiqc of alignment results
    SamStatsQC(ch_indexedBam, inName)
}
