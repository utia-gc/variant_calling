/**
 * Process to run GATK GenomicsDBImport
 *
 * Import single-sample GVCFs into GenomicsDB before joint genotyping.
 * @see https://gatk.broadinstitute.org/hc/en-us/articles/13832686645787-GenomicsDBImport
 *
 * @input gvcfs a channel of a collection of GVCF files.
 * @input reference the reference sequence file channel of format [fasta, fai, dict].
 * @emit genomicsDB the GenomicsDB of merged GVCF files and interval used of foramt [interval, GenomicsDB].
 */
process gatk_GenomicsDBImport {
    label 'gatk'

    input:
        path gvcfs
        path gvcfTbis
        tuple path(fasta), path(fai), path(dict)
        val interval

    output:
        tuple val(interval), path('*gdb'), emit: genomicsDB

    script:
        // create the string of INPUT arguments
        // prepend '--INPUT' to each BAM
        // join into a space-delimted string
        def variants = gvcfs.collect { gvcf ->
            "--variant ${gvcf}" 
        }.join(' ')

        String args = new Args(task.ext).buildArgsString()

        String javaOptions = ""

        """
        mkdir \${PWD}/tmp

        gatk --java-options "${javaOptions}" GenomicsDBImport \\
            ${variants} \\
            --genomicsdb-workspace-path ${interval}_gdb \\
            --reference ${fasta} \\
            --intervals ${interval} \\
            ${args}
        """
}
