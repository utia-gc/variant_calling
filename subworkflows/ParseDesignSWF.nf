include { ParseDesign   } from "${projectDir}/modules/ParseDesign.nf"

workflow ParseDesignSWF {
    take:
        design

    main:
        ParseDesign(
            design
        )
            .csv
            .splitCsv(header:true, sep:',')
            .map { createDesignChannel(it) }
            .branch {
                reads: it[1].any { it =~ /(fastq|fq)/ }   // add channels with fastq files (either fastq or fq) to reads channel
            }
            .set { design }

    emit:
        reads  = design.reads
}


// create a list of data from the csv
def createDesignChannel(LinkedHashMap row) {
    // reads
    if (row.fq1) {
        // store metadata in a Map
        def metadata = [:]
        metadata.libID      = row.lib_ID
        metadata.sampleName = row.sample_rep
        metadata.readType   = row.fq2 ? 'paired' : 'single'

        // store reads in a list
        def reads = []
        if (metadata.readType == 'single') {
            reads = [file(row.fq1)]
        } else {
            reads = [file(row.fq1), file(row.fq2)]
        }

        return [metadata, reads]
    }
}
