/**
 * Build the stem name for a sample (e.g. reads or alignments) from a map of metadata.
 *
 * This is just a way to reuse code that dynamically adds sample number and lane info to sample names if they exist.
 *
 * @param metadata A LinkedHashMap of metadata for a sample.
 */
static String buildStemName(LinkedHashMap metadata) {
    // the stem name must start with the sample name
    ArrayList stemNameComponents = [metadata.sampleName]

    // add lane number to stem name if it exists
    if(metadata.lane) stemNameComponents += "L${metadata.lane}"

    return stemNameComponents.join('_')
}
