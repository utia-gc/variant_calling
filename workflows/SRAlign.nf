/*
=====================================================================
    SRAlign WORKFLOW
=====================================================================
*/

/*
This object takes care of many necessary steps upon construction:
    - Logs a header for the pipeline that prints pipeline name and logo
    - Prints a help message if help parameter is specified
    - Checks parameters
*/ 
def srawf = new SRAlignWorkflow(log, params, workflow)


/*
    ---------------------------------------------------------------------
    Set useful pipeline values
    ---------------------------------------------------------------------
*/

// set output filename base prefix
outBasePrefix   = srawf.outBasePrefix
outUniquePrefix = srawf.outUniquePrefix

// set genome and contaminant values
genome = params.genomes[ params.genome ]
contaminant = params.genomes[ params.contaminant ]

/*
    ---------------------------------------------------------------------
    Import modules
    ---------------------------------------------------------------------
*/

include { ParseDesignSWF        as ParseDesign        } from "${projectDir}/subworkflows/inputs/ParseDesignSWF.nf"
include { TrimReadsSWF          as TrimReads          } from "${baseDir}/subworkflows/reads/TrimReadsSWF.nf"
include { ReadsQCSWF            as ReadsQC            } from "${baseDir}/subworkflows/reads/ReadsQCSWF.nf"
include { AlignBowtie2SWF       as AlignBowtie2       ; 
          AlignBowtie2SWF       as ContamBowtie2      } from "${baseDir}/subworkflows/align/AlignBowtie2SWF.nf"
include { AlignHisat2SWF        as AlignHisat2        ; 
          AlignHisat2SWF        as ContamHisat2       } from "${baseDir}/subworkflows/align/AlignHisat2SWF.nf"
include { PreprocessSamSWF      as PreprocessSam      } from "${baseDir}/subworkflows/align/PreprocessSamSWF.nf"
include { SamStatsQCSWF         as SamStatsQC         } from "${baseDir}/subworkflows/align/SamStatsQCSWF.nf"
include { SeqtkSample           as SeqtkSample        } from "${baseDir}/modules/reads/SeqtkSample.nf"
include { ContaminantStatsQCSWF as ContaminantStatsQC } from "${baseDir}/subworkflows/align/ContaminantStatsQCSWF.nf"
include { PreseqSWF             as Preseq             } from "${baseDir}/subworkflows/align/PreseqSWF.nf"
include { DeepToolsMultiBamSWF  as DeepToolsMultiBam  } from "${projectDir}/subworkflows/align/DeepToolsMultiBamSWF.nf"
include { FullMultiQC           as FullMultiQC        } from "${baseDir}/modules/misc/FullMultiQC.nf"


workflow SRAlign {
    /*
    ---------------------------------------------------------------------
        Read design file, parse sample names and identifiers, and stage reads files
    ---------------------------------------------------------------------
    */

    // set channel for input design file
    ch_input = file(params.input)

    // Subworkflow: Parse design file
    ParseDesign(
        ch_input
    )
    ch_rawReads         = ParseDesign.out.reads
    ch_bamIndexedGenome = ParseDesign.out.bamBai


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
    } else {
        ch_trimReads = Channel.empty()
    }


    /*
    ---------------------------------------------------------------------
        Reads QC
    ---------------------------------------------------------------------
    */

    if (!params.skipReadsQC) {
        // Subworkflow: FastQC and MulitQC for raw and trimmed reads
        ReadsQC(
            ch_rawReads,
            ch_trimReads,
            outUniquePrefix
        )
        ch_rawReadsFQC  = ReadsQC.out.raw_fqc_zip
        ch_trimReadsFQC = ReadsQC.out.trim_fqc_zip
    } else {
        ch_rawReadsFQC  = Channel.empty()
        ch_trimReadsFQC = Channel.empty()
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
        ch_samGenome = Channel.empty()
        // Align reads to genome
        switch (params.alignmentTool) {
            case 'bowtie2':
                // Subworkflow: Align reads to genome with bowtie2 and build index if necessary
                AlignBowtie2(
                    ch_readsToAlign,
                    genome,
                    params.genome
                )
                ch_samGenome = AlignBowtie2.out.sam
                break
            
            case 'hisat2':
                // Subworkflow: Align reads to genome with hisat2 and build index if necessary
                AlignHisat2(
                    ch_readsToAlign,
                    genome,
                    params.genome,
                    params.forceUseHisat2Index,
                    params.buildSpliceAwareIndex
                )
                ch_samGenome = AlignHisat2.out.sam
                break
    }

    // Preprocess sam files: mark duplicates, sort alignments, compress to bam, and index
    PreprocessSam(
        ch_samGenome
    )
    ch_bamIndexedGenome = PreprocessSam.out.bamBai.mix(ch_bamIndexedGenome)


    if (!params.skipSamStatsQC) {
        // Subworkflow: Samtools stats and samtools idxstats and multiqc of alignment results
        SamStatsQC(
            ch_bamIndexedGenome,
            outUniquePrefix
        )
        ch_alignGenomeStats    = SamStatsQC.out.samtoolsStats
        ch_alignGenomeIdxstats = SamStatsQC.out.samtoolsIdxstats
        ch_alignGenomePctDup   = SamStatsQC.out.pctDup
        }
    } else {
        ch_alignGenomeStats    = Channel.empty()
        ch_alignGenomeIdxstats = Channel.empty()
        ch_alignGenomePctDup   = Channel.empty()
    }

    /*
    ---------------------------------------------------------------------
        Check contamination 
    ---------------------------------------------------------------------
    */

    if (params.contaminant && !params.skipAlignContam) {
        ch_samContaminant = Channel.empty()

        // Sample reads
        SeqtkSample(
            ch_readsToAlign
        ) 
        ch_readsContaminant = SeqtkSample.out.sampleReads

        // Align reads to contaminant genome
        switch (params.alignmentTool) {
            case 'bowtie2':
                // Subworkflow: Align reads to contaminant genome with bowtie2 and build index if necessary
                ContamBowtie2(
                    ch_readsContaminant,
                    contaminant,
                    params.contaminant
                )
                ch_samContaminant = ContamBowtie2.out.sam
                break
            
            case 'hisat2':
                // Subworkflow: Align reads to contaminant genome with hisat2 and build index if necessary
                ContamHisat2(
                    ch_readsContaminant,
                    contaminant,
                    params.contaminant,
                    params.forceUseHisat2Index,
                    false
                )
                ch_samContaminant = ContamHisat2.out.sam
                break
        }

        // Get contaminant alignment stats
        ContaminantStatsQC(
            ch_samContaminant,
            outUniquePrefix
        )
        ch_contaminantFlagstat = ContaminantStatsQC.out.samtoolsFlagstat
    } else {
        ch_contaminantFlagstat = Channel.empty()
    }

    /*
    ---------------------------------------------------------------------
        Dataset stats
    ---------------------------------------------------------------------
    */

    // Preseq
    if (!params.skipPreseq) {
        Preseq(
            ch_bamIndexedGenome
        )
        ch_preseqLcExtrap = Preseq.out.psL
    } else {
        ch_preseqLcExtrap = Channel.empty()
    }

    // deepTools
    ch_alignments = ch_bamIndexedGenome
    ch_alignmentsCollect = 
        ch_alignments
        .multiMap {
            it ->
            bam:     it[1]
            bai:     it[2]
            toolIDs: it[3]
        }

    DeepToolsMultiBam(
        ch_alignmentsCollect.bam.collect(),
        ch_alignmentsCollect.bai.collect(),
        outBasePrefix
    )
    ch_corMatrix = DeepToolsMultiBam.out.corMatrix
    ch_PCAMatrix = DeepToolsMultiBam.out.PCAMatrix


    /*
    ---------------------------------------------------------------------
        Full pipeline MultiQC
    ---------------------------------------------------------------------
    */

    ch_fullMultiQC = Channel.empty()
        .concat(ch_rawReadsFQC)
        .concat(ch_trimReadsFQC)
        .concat(ch_alignGenomeStats)
        .concat(ch_alignGenomeIdxstats)
        .concat(ch_alignGenomePctDup)
        .concat(ch_contaminantFlagstat)
        .concat(ch_preseqLcExtrap)
        .concat(ch_corMatrix)
        .concat(ch_PCAMatrix)
    
    // set channel for MultiQC config file
    ch_multiqcConfig = file(params.multiqcConfig)

    FullMultiQC(
        outUniquePrefix,
        ch_multiqcConfig,
        ch_fullMultiQC.collect()
    )
}
