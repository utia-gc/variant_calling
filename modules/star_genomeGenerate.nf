process star_genomeGenerate {
    tag "${genome.name}"

    label 'star'

    label 'big_cpu'
    label 'sup_mem'
    label 'big_time'

    input:
        path genome
        path annotationsGTF

    output:
        path 'STAR', emit: index

    script:
        String args = new Args(task.ext).buildArgsString()

        """
        mkdir STAR

        STAR \
            --runMode genomeGenerate \
            --genomeFastaFiles ${genome} \
            --sjdbGTFfile ${annotationsGTF} \
            --genomeDir STAR \
            --runThreadN ${task.cpus} \
            ${args}
        """
}
