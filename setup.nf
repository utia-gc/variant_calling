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

workflow {
    readsSource = file(params.readsSource)
    readsDest   = file(params.readsDest)

    COPY_READS(
        readsSource,
        readsDest
    )
}
