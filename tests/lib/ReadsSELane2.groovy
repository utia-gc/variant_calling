import groovy.transform.InheritConstructors

@InheritConstructors
class ReadsSELane2 extends Reads {
    LinkedHashMap metadata = [
        sampleName:   'SRR1066657',
        sampleNumber: '3',
        lane:         '002',
        readType:     'single'
    ]
    List reads = [
        'https://github.com/utia-gc/ngs-test/raw/ngs/data/reads/raw/SRR1066657_S3_L002_R1_001.fastq.gz',
        'https://github.com/utia-gc/ngs-test/raw/ngs/data/reads/raw/SRR1066657_S3_L002_R1_001.fastq.gz.NOFILE'
    ]
}
