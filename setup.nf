workflow {
    /*
    * Work with directory that currently has reads
    */
    println "Reads directory: ${params.readsDir}"
    def readsDir = file(params.readsDir)

    // define list of patterns to match
    def patterns = [
        /.*/
    ]
    // combine patterns into single regex
    def combinedPattern = patterns.join('|')

    /*
    * Work with directory to copy reads to
    */
    println "Copy directory: ${params.copyDir}"
    def copyDir = file(params.copyDir)

    // make the copy dir
    copyDir.mkdirs()

    // iterate through files in directory that match combined regex
    readsDir.eachFileMatch(~combinedPattern) { fastq ->
        println "Matching File: ${fastq}"
        // copy files to copy dir
        fastq.copyTo(copyDir)
    }
}
