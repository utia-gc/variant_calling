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
