abstract class Reads {
    public getSampleName() {
        this.metadata.sampleName
    }
    
    public getSampleNumber() {
        this.metadata.sampleNumber
    }
    
    public getLane() {
        this.metadata.lane
    }
    
    public getReadType() {
        this.metadata.readType
    }
    
    public getR1() {
        this.reads[0]
    }
    
    public getR2() {
        this.reads[1]
    }

    /**
     * Get the stem name of a fastq file.
     * 
     * Get the stem name of a fastq file. 
     * This is the part of the fastq file preceding '_R1_001'.
     * For fastq files that follow Illumina naming conventions this should be the same as '<SampleName>_S<SampleNumber>_L<LaneNumber>'.
     * @see https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
     *
     * @return String fastq stem name.
     */
    public getStemName() {
        file(this.getR1()).getName() - ~/_R1_001.*/
    }

    public getR1SimpleName() {
        // regexp pattern that matches common fastq filename endings
        // matches: fastq.gz, fq.gz, fastq, fq
        def fastqSuffix = ~/\.f(?:ast)?q(?:\.gz)?$/

        file(this.getR1()).getName() - fastqSuffix
    }

    public getR2SimpleName() {
        // regexp pattern that matches common fastq filename endings
        // matches: fastq.gz, fq.gz, fastq, fq
        def fastqSuffix = ~/\.f(?:ast)?q(?:\.gz)?$/

        file(this.getR2()).getName() - fastqSuffix
    }
}
