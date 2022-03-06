include { Bowtie2Build } from '../modules/Bowtie2Build.nf'
include { AlignBowtie2 } from '../modules/AlignBowtie2.nf'

workflow AlignBowtie2SWF {
    take:
        reads
    
    main:
        bowtie2Index = params.genomes[params.genome].bowtie2 ?: false
        if (bowtie2Index) {
            bowtie2Indexes = Channel
                .fromPath("${bowtie2Index}*", checkIfExists: true)
                .collect()
        } else {
            Bowtie2Build(
                params.genomes[params.genome].fasta
            )
            .collect()
            .set { bowtie2Indexes }
        }

        AlignBowtie2(reads, bowtie2Indexes)

    emit:
        sam = AlignBowtie2.out.sam
}
