workflow {
    readsSource = file(params.readsSource)
    readsDest   = file(params.readsDest)

    COPY_READS(
        readsSource,
        readsDest
    )

    samplesheet = file(params.samplesheet)

    WRITE_SAMPLESHEET(
        readsDest,
        samplesheet
    )
}

workflow COPY_READS {
    take:
        source_reads_dir
        destination_reads_dir

    main:
        // make the destination directory
        destination_reads_dir.mkdirs()

        // iterate through all fastq.gz files in source directory
        source_reads_dir.eachFileMatch(~/.*\.fastq\.gz/) { fastq ->
            // copy files to copy dir
            fastq.copyTo(destination_reads_dir)
        }
}

workflow WRITE_SAMPLESHEET {
    take:
        reads_dir
        samplesheet

    main:
        // create channel of read pairs
        Channel
            .fromFilePairs("${reads_dir}/*_R{1,2}_001.fastq.gz", size: -1, checkIfExists: true)
            .set { ch_readPairs }
        ch_readPairs.dump(tag: "ch_readPairs", pretty: true)

        // cast ch_readPairs to a map and write to a file
        ch_readPairs
            .map { stemName, reads ->
                "${stemName},${reads[0]},${reads[1]}"
            }
            .collectFile(
                name: samplesheet.name,
                newLine: true,
                storeDir: samplesheet.parent,
                sort: true,
                seed: 'sampleName,reads1,reads2'
            )
}
