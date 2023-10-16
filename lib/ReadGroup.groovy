import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

/**
 * Build the read group line from sample metadata and matcher data from the sequence identifier.
 *
 * @params LinkedHashMap metadata A metadata map.
 * @params Matcher matcher A matcher object of sequence identifier.
 *
 * @return String A tab-separated formatted RG line
 */
static String buildRGLine(metadata, matcher) {
    ArrayList rgFields = ['@RG']
    def rgID = ''
    def rgPU = ''

    // add dynamically determined read group fields
    // this includes ID and PU
    if(matcher.find()) {
        rgID = "${metadata.sampleName}_${matcher.group('instrument')}_${matcher.group('runNumber')}_${matcher.group('flowcellID')}.${matcher.group('lane')}"
        rgPU = rgID
    } else {
        rgID = "${metadata.sampleName}.${metadata.lane}"
    }
    rgFields += "ID:${rgID}"
    rgPU ? rgFields += "PU:${rgPU}" : null

    // add more straightforwardly determined fields
    rgFields += "SM:${metadata.sampleName}"
    rgFields += "LB:${metadata.sampleName}"

    // add static fields
    rgFields += "PL:ILLUMINA"

    return rgFields.join('\t')
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
