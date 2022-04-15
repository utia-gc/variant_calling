include { ParseDesign } from "${projectDir}/modules/inputs/ParseDesign.nf"

workflow ParseDesignBamsSWF{
    take:
        design

    main:
        ParseDesign(
            'alignments',
            design
        )
            .csv
            .splitCsv(header:true, sep:',')
            .map { createBamsChannel(it) }
            .set { bams }

    emit:
        bams 
}


// create a list of data from the csv
def createBamsChannel(LinkedHashMap row) {
    // store metadata in a Map
    def metadata = [:]
    metadata.libID      = row.lib_ID
    metadata.sampleName = row.sample_rep

    // store bams
    bam = row.bam

    // create an empty list for tool IDs for suffixes
    toolIDs = row.tool_IDs.split('_')

    // check that reads files exist
    bam.each {
        if (!it.exists()) {
            exit 1, "ERROR: ${it} does not exist!"
        }
    }

    return [metadata, bam, toolIDs]
}
