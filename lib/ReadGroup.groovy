import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

static def buildBwaMem2RGLine(rgFields) {
    ArrayList rgLineElements = ['@RG']

    rgFields.each { tag, value ->
        rgLineElements += "${tag}:${value}"
    }

    return rgLineElements.join('\t')
}

/**
 * Build the read group fields from sample metadata and matcher data from the sequence identifier.
 *
 * @params LinkedHashMap metadata A metadata map.
 * @params Matcher matcher A matcher object of sequence identifier.
 *
 * @return LinkedHashMap Read group fields containing at least ID, SM, LB, and PL fields.
 */
static LinkedHashMap buildRGFields(metadata, matcher) {
    def rgFields = [:]

    // add dynamically determined read group fields -- ID and PU
    if(matcher.find()) {
        rgFields += ['ID': "${matcher.group('instrument')}_${matcher.group('runNumber')}_${matcher.group('flowcellID')}.${matcher.group('lane')}"]
        rgFields += ['PU': "${matcher.group('flowcellID')}.${matcher.group('lane')}.${matcher.group('index')}"]
    } else {
        rgFields += ['ID': "${metadata.sampleName}.${metadata.lane}"]
    }

    // add more straightforwardly determined dynamic fields
    rgFields += ['SM': "${metadata.sampleName}"]
    rgFields += ['LB': "${metadata.sampleName}"]

    // add static fields
    rgFields += ['PL': "ILLUMINA"]

    return rgFields
}

/**
 * Match sequence identifier against valid Illumina fastq file format
 *
 * @params String sequenceIdentifier The sequence identifier line. Should have the leading '@' stripped.
 *
 * @return Matcher A matcher for the sequenceIdentifier against a validated sequence identifer regex pattern.
 */
static def matchSequenceIdentifier(sequenceIdentifier) {
    def validSequenceIdentifierPattern = /^(?<instrument>[a-zA-Z0-9_]+):(?<runNumber>[0-9]+):(?<flowcellID>[a-zA-Z0-9]+):(?<lane>[0-9]+):(?<tile>[0-9]+):(?<xPos>[0-9:]+):(?<yPos>[0-9]+)(:(?<umi>[ATGCN+]+))? (?<read>[12]):(?<isFiltered>[YN]):(?<controlNumber>[0-9]+):(?<index>[ATGCN+]+)$/

    return (sequenceIdentifier =~ validSequenceIdentifierPattern)
}

/**
 * Read the first sequence identifier from a fastq.gz file
 *
 * @params fastqPath Path object to fastq file.
 *
 * @return String sequence identifier line.
 */
static String readFastqFirstSequenceIdentifier(fastqPath) {
    def sequenceIdentifier = ''
    def fastqInputStream = new GZIPInputStream(fastqPath.newInputStream())
    try {
        def reader = new BufferedReader(new InputStreamReader(fastqInputStream, 'UTF-8'))
        sequenceIdentifier = reader.readLine()[1..-1]
    } finally {
        fastqInputStream.close()
    }

    return sequenceIdentifier
}
