/*
---------------------------------------------------------------------
    HANDLE INPUTS
---------------------------------------------------------------------
*/

/*
    Design and Inputs
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
    Genomes and References
*/

// check genome
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "Genome '${params.genome}' is not available.\n\tAvailable genomes include: ${params.genomes.keySet().join(", ")}"
}


// set genome variables
bt2Index = params.genome ? params.genomes[params.genome].bowtie2 ?: false : false


/*
---------------------------------------------------------------------
    MAIN WORKFLOW
---------------------------------------------------------------------
*/

include { ParseDesignSWF  as ParseDesign  } from '../subworkflows/ParseDesignSWF.nf'
include { RawReadsQCSWF   as RawReadsQC   } from '../subworkflows/RawReadsQCSWF.nf'
include { TrimReadsSWF    as TrimReads    } from '../subworkflows/TrimReadsSWF.nf'
include { TrimReadsQCSWF  as TrimReadsQC  } from '../subworkflows/TrimReadsQCSWF.nf'
include { AlignBowtie2SWF as AlignBowtie2 } from '../subworkflows/AlignBowtie2SWF.nf'
include { SamStatsQCSWF   as SamStatsQC   } from '../subworkflows/SamStatsQCSWF.nf'


workflow sralign {
    ParseDesign(ch_input)
    RawReadsQC(ParseDesign.out.rawReads, inName)
    TrimReads(ParseDesign.out.rawReads)
    TrimReadsQC(TrimReads.out.trimReads, inName) 

    AlignBowtie2(TrimReads.out.trimReads, bt2Index)
    SamStatsQC(AlignBowtie2.out.bamBai, inName)

}
