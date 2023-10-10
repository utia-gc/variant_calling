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
    // fill in required metadata
    metadata.sampleName = row.sampleName
    metadata.readType   = row.reads2 ? "paired" : "single"
    // add optional metadata
    if (row.sampleNumber) {
        metadata.sampleNumber = row.sampleNumber
    }
    if (row.lane) {
        metadata.lane = row.lane
    }

    // store reads in lists
    def reads1 = file(row.reads1)
    if(metadata.readType == "single") {
        def emptyFileName = "${reads1.name}.NOFILE"
        file("${projectDir}/assets/NO_FILE").copyTo(emptyFileName)
        reads2 = file(emptyFileName)
    } else if(metadata.readType == 'paired') {
        reads2 = file(row.reads2)
    }

    return [metadata, [reads1], [reads2]]
}
