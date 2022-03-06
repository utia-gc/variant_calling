/*
=====================================================================
    HANDLE INPUTS
=====================================================================
*/


/*
    ---------------------------------------------------------------------
    Tools
    ---------------------------------------------------------------------
*/

def tools = [
    trim      : ['fastp'],
    alignment : ['bowtie2', 'hisat2']
]

// check valid read-trimming tool
assert params.trimTool in tools.trim , 
    "'${params.trimTool}' is not a valid read trimming tool.\n\tValid options: ${tools.trim.join(', ')}\n\t" 

// check valid alignment tool
assert params.alignmentTool in tools.alignment , 
    "'${params.alignmentTool}' is not a valid alignment tool.\n\tValid options: ${tools.alignment.join(', ')}\n\t"


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


/*
=====================================================================
    MAIN WORKFLOW
=====================================================================
*/

include { ParseDesignSWF   as ParseDesign   } from '../subworkflows/ParseDesignSWF.nf'
include { RawReadsQCSWF    as RawReadsQC    } from '../subworkflows/RawReadsQCSWF.nf'
include { TrimReadsSWF     as TrimReads     } from '../subworkflows/TrimReadsSWF.nf'
include { TrimReadsQCSWF   as TrimReadsQC   } from '../subworkflows/TrimReadsQCSWF.nf'
include { AlignBowtie2SWF  as AlignBowtie2  } from '../subworkflows/AlignBowtie2SWF.nf'
include { PreprocessSamSWF as PreprocessSam } from '../subworkflows/PreprocessSamSWF.nf'
include { SamStatsQCSWF    as SamStatsQC    } from '../subworkflows/SamStatsQCSWF.nf'


workflow sralign {

    /*
    ---------------------------------------------------------------------
        Read design file, parse sample names and identifiers, and stage reads files
    ---------------------------------------------------------------------
    */

    // Subworkflow: Parse design file
    ParseDesign(
        ch_input
    )
    ch_rawReads = ParseDesign.out.rawReads


    /*
    ---------------------------------------------------------------------
        Raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipRawFastQC) {
        // Subworkflow: Raw reads fastqc and mulitqc
        RawReadsQC(
            ch_rawReads,
            inName
        )
    }


    /*
    ---------------------------------------------------------------------
        Trim raw reads
    ---------------------------------------------------------------------
    */

    if (!params.skipTrimReads) {
        // Trim reads
        switch (params.trimTool) {
            case 'fastp':
                // Subworkflow: Trim raw reads
                TrimReads(
                    ch_rawReads
                )
                ch_trimReads = TrimReads.out.trimReads
                break
        }


        // Trimmed reads QC
        if (!params.skipTrimReadsQC) {
            // Subworkflow: Trimmed reads fastqc and multiqc
            TrimReadsQC(
                ch_trimReads,
                inName
            ) 
        }
    } 


    /*
    ---------------------------------------------------------------------
        Align reads to genome
    ---------------------------------------------------------------------
    */

    // Set channel of reads to align 
    if (!params.forceAlignRawReads) {
        if (!params.skipTrimReads) {
            ch_readsToAlign = ch_trimReads
        } else {
            ch_readsToAlign = ch_rawReads
        }
    } else {
        ch_readsToAlign = ch_rawReads
    }


    if (!params.skipAlignGenome) {
        // Align reads to genome
        switch (params.alignmentTool) {
            case 'bowtie2':
                // Subworkflow: Align reads to genome with bowtie2 and build index if necessary
                AlignBowtie2(
                    ch_readsToAlign
                )
                ch_samGenome = AlignBowtie2.out.sam
                break
        }
    

    PreprocessSam(
        ch_samGenome
    )
    ch_bamGenome        = PreprocessSam.out.bam
    ch_bamIndexedGenome = PreprocessSam.out.bamBai


        if (!params.skipSamStatsQC) {
            // Subworkflow: Samtools stats and samtools idxstats and multiqc of alignment results
            SamStatsQC(
                ch_bamIndexedGenome,
                inName
            )
        }
    }
}
