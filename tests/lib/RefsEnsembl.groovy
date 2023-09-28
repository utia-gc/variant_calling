import groovy.transform.InheritConstructors

@InheritConstructors
class RefsEnsembl extends Refs {
    String name = 'R64-1-1'
    LinkedHashMap sequences = [
        genome: [
            decomp: 'https://github.com/utia-gc/ngs-test/raw/ngs/data/references/R64-1-1/genome_I.fa',
            gzip:   'https://github.com/utia-gc/ngs-test/raw/ngs/data/references/R64-1-1/genome_I.fa.gz',
            index:  'https://github.com/utia-gc/ngs-test/raw/ngs/data/references/R64-1-1/genome_I.fa.fai'
        ]
    ]
    LinkedHashMap annotations = [
        gtf: [
            decomp: 'https://github.com/utia-gc/ngs-test/raw/ngs/data/references/R64-1-1/annotations_I.gtf',
            gzip:   'https://github.com/utia-gc/ngs-test/raw/ngs/data/references/R64-1-1/annotations_I.gtf.gz'
        ]
    ]
}
