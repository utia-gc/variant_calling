class SRAlignWorkflow {
    /*
    ---------------------------------------------------------------------
        Fields and properties
    ---------------------------------------------------------------------
    */

    // workflow features
    public static def log

    // workflow name
    public static final String pipelineName = "trev-f/SRAlign"

    // workflow title in ASCII characters
    public static final String pipelineTitle = (
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


    /*
    ---------------------------------------------------------------------
        Methods
    ---------------------------------------------------------------------
    */
        
    // constructor method
    SRAlignWorkflow(log) {
        // set log
        this.log = log

        // display the header
        log.info createHeader()
    }

    // create a header
    public static String createHeader() {
        // initialize an empty header list
        def header = []

        // add pipeline header info to list
        header.add(pipelineTitle)
        header.add(pipelineName)

        // return header info as string with new line breaks
        header.join("\n")
    }
}
