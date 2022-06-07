import groovy.json.JsonSlurper
/**
 * SRAlignWorkflow class is a class to perform necessary workflow functions upon initialization.
 * These functions include: creating a header for STDOUT and documentation, generating help documentation, checking that specified params are valid, and others.
*/
class SRAlignWorkflow {
    /*
    ---------------------------------------------------------------------
        Fields and properties
    ---------------------------------------------------------------------
    */

    // workflow features
    /** built-in Nextflow logging */
    public static def           log
    /** Parameters specified in the workflow */
    public static def           params
    /** Specifications for parameters such as default values, param descriptions, etc. to be used in creating help documentation and check parameters. */
    public static LinkedHashMap paramSpecs

    /** valid tools available in the pipeline */
    public static LinkedHashMap validTools = [
        trim      : ['fastp'],
        alignment : ['bowtie2', 'hisat2']
    ]

    // pipeline ASCII logo
    public static final String pipelineLogo = (
        """
         ######  ########     ###    ##       ####  ######   ##    ## 
        ##    ## ##     ##   ## ##   ##        ##  ##    ##  ###   ## 
        ##       ##     ##  ##   ##  ##        ##  ##        ####  ## 
         ######  ########  ##     ## ##        ##  ##   #### ## ## ## 
              ## ##   ##   ######### ##        ##  ##    ##  ##  #### 
        ##    ## ##    ##  ##     ## ##        ##  ##    ##  ##   ### 
         ######  ##     ## ##     ## ######## ####  ######   ##    ## 
        """
    )

    // pipeline name
    public static final String pipelineName = "trev-f/SRAlign"

    // purpose statement
    public static final String purpose = (
        """
        A flexible pipeline for short read alignment to a reference with extensive QC reporting.
        """
    )

    // general nextflow help
    public static final String generalNFHelp = (
        """
        Nextflow command-line options:
            nextflow help

        Nextflow run command-line options:
            nextflow run -help
        """
    )

    // usage statement
    public static final String usage = (
        """
        nextflow run trev-f/SRAlign -profile docker --input <input.csv> --genome <valid genome key>
        """
    )

    /*
    ---------------------------------------------------------------------
        Methods
    ---------------------------------------------------------------------
    */
        
    /**
     * SRAlignWorkflow constructor method
     *
     * @param log the Nextflow log object
     * @param params the workflow parameters
     * @param paramSpecs the parameters specifications and defaults
     *
     * @return SRAlignWorkflow object
    */
    SRAlignWorkflow(log, params) {
        // set log
        this.log        = log
        this.params     = params
        this.paramSpecs = (new JsonSlurper()).parse(new File('/home/treevooor/SRAlign/parameter_specifications.json'))

        // add options to paramSpecs
        LinkedHashMap paramSpecs = addValidOptions(paramSpecs)

        // display the header
        log.info createHeader()

        // display a help message and exit the program
        if (params.help) {
            log.info getHelp(paramSpecs)
            System.exit(0)
        }
    }


    /**
     * Creates a header from pipeline info to be displayed on STDOUT and in logging.
     *
     * @return a header message
    */
    public static String createHeader() {
        // initialize an empty header list
        def header = []

        // add pipeline header info to list
        header.add(pipelineLogo)
        header.add(pipelineName)
        header.add(purpose)
        header.add("")              // extra newline at end helps formatting look better

        // return header info as string with new line breaks
        header.join("\n")
    }


    /**
     * Creates a method for adding options to paramSpecs
     *
     * @param paramSpecs parameter specifications
     *
     * @return parameter specifications with valid options added
    */
    public static LinkedHashMap addValidOptions(paramSpecs) {
        // add valid genome options to genome and contaminant
        paramSpecs.genome.options        = params.genomes.keySet()
        paramSpecs.contaminant.options   = params.genomes.keySet()

        // add valid tools
        paramSpecs.trimTool.options      = validTools.trim
        paramSpecs.alignmentTool.options = validTools.alignment

        return paramSpecs
    }


    /**
     * Creates a help string for a set of parameter specifications.
     * This is a helper function to be used with `writeHelpDoc` to write help documentation.
     *
     * @param specs a set of parameter specifications
     *
     * @return a help message string for a set of parameter specifications
    */
    public static String writeHelpString(LinkedHashMap specs) {
        // intialize an empty documentation list
        List doc = []
        
        // iterate through parameters
        specs.each {
            spec ->
            // add parameter name and description to help doc as a single line
            doc.add("--${spec.key.padRight(30)}${spec.value.description}")

            // add parameter default value to help doc
            if (spec.value.default) {
                doc.add("${"".padRight(32)}Default: ${spec.value.default}")
            }

            // add parameter options to help doc
            if (spec.value.options) {
                doc.add("${"".padRight(32)}Options: ${spec.value.options.join(", ")}")
            }

            // extra newline at end helps formatting look better
            doc.add("")
        }

        // join help docs together into a string all on new lines
        doc.join("\n")
    }


    /**
     * Creates a mature string of help documentation for parameter specifications.
     * This method splits up parameter help docs into three sets: required params, optional params, and params that skip steps of the workflow.
     *
     * @param specs parameter specifications
     *
     * @return a mature help message string for parameters
    */
    public static String writeParamHelpDoc(LinkedHashMap specs) {
        // initialize an emtpy list of help docs
        List help = []

        // add required params
        help.add("Required parameters:")
        help.add(
            writeHelpString(
                specs.findAll {
                    it.value.required == true
                }
            )
        )

        // add optional params
        help.add("Optional parameters:")
        help.add(
            writeHelpString(
                specs.findAll {
                    it.value.required != true && it.value.skip != true
                }
            )
        )

        // add skip params
        help.add("Skip module parameters:")
        help.add(
            writeHelpString(
                specs.findAll {
                    it.value.required != true && it.value.skip == true
                }
            )
        )

        // extra newline helps formatting look better
        help.join("\n")
    }


    /**
     * Creates a help message for parameter specifications
     *
     * @param paramSpecs the parameters specifications and defaults
     *
     * @return help documentation
    */
    public static String getHelp(paramSpecs) {
        // initialize an empty help list
        def help = []

        // add help message info to list
        help.add("Nextflow command-line help:${generalNFHelp}")
        help.add("")                // extra newline helps formatting look better
        help.add("Usage:${usage}")
        help.add("")
        help.add(writeParamHelpDoc(paramSpecs))
        help.add("")

        // return help message as string with new line breaks
        help.join("\n")
    }
}
