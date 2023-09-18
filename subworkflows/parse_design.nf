workflow Parse_Design {
    take:
        samplesheet

    main:
        Channel
            .fromPath( samplesheet, checkIfExists: true )
            .splitCsv( header:true, sep:',' )
            .map { createSampleReadsChannel(it) }
            .set { ch_samples }

    emit:
        samples = ch_samples
}


// create a list of data from the csv
def createSampleReadsChannel(LinkedHashMap row) {
    // store metadata in a Map
    def metadata = [:]
    metadata.sampleName = row.sampleName
    metadata.readType   = row.reads2 ? "paired" : "single"

    // store reads in a list
    def reads = []
    if(metadata.readType == "single") {
        reads = [file(row.reads1)]
    } else {
        reads = [file(row.reads1), file(row.reads2)]
    }

    return [metadata, reads]
}