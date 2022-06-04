class SRAlignWorkflow {
    /*
    ---------------------------------------------------------------------
        Fields and properties
    ---------------------------------------------------------------------
    */

    // workflow features
    public static def log
    public static def params

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
        
    // constructor method
    SRAlignWorkflow(log, params) {
        // set log
        this.log    = log
        this.params = params

        // display the header
        log.info createHeader()

        // display a help message and exit the program
        if (params.help) {
            log.info getHelp()
            System.exit(0)
        }
    }


    // create a header
    public static String createHeader() {
        // initialize an empty header list
        def header = []

        // add pipeline header info to list
        header.add(pipelineLogo)
        header.add(pipelineName)
        header.add("")            // extra newline at end helps formatting look better

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
