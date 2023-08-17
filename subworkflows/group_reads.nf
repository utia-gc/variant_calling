workflow Group_Reads {
    take:
        reads

    main:
        reads
            .map { selectReadsMetadataMerging(it) }
            .groupTuple()
            .map { reshapeGroupedLaneReads(it) }
            .dump(tag: "ch_reads_grouped", pretty: true)
            .set { ch_reads_grouped }

    emit:
        reads_grouped = ch_reads_grouped
}

def selectReadsMetadataMerging(ArrayList reads) {
    def metadata = [:]
    metadata.put("sampleName", reads[0].get("sampleName"))
    metadata.put("readType", reads[0].get("readType"))

    reads[0] = metadata

    return reads
}

def reshapeGroupedLaneReads(groupedReads) {
    // make empty lists for R1 and R2
    def reads1 = []
    def reads2 = []

    groupedReads[1].eachWithIndex { it, index ->
        reads1.add(it[0])
        reads2.add(it[1] ?: "${projectDir}/assets/NO_FILE${index}")
    }

    return [groupedReads[0], reads1.sort(), reads2.sort()]
}
