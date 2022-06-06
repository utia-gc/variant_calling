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

    // usage statement
    public static final String usage = (
        """
        nextflow run trev-f/SRAlign -profile docker --input input.csv --genome WBcel235
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
    SRAlignWorkflow(log, params, paramSpecs) {
        // set log
        this.log        = log
        this.params     = params
        this.paramSpecs = paramSpecs

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


    // create a help message
    public static String getHelp() {
        // initialize an empty help list
        def help = []

        // add help message info to list
        help.add("Usage:${usage}")
        help.add("")                // extra newline at end helps formatting look better

        // return help message as string with new line breaks
        help.join("\n")
    }
}
