class SRAlignWorkflow {
    // import log
    public static def          log

    // workflow name
    public static final String workflowName = "trev-f/SRAlign"

    // workflow title in ASCII characters
    public static final String workflowHeader = (
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
        
    // constructor method
    SRAlignWorkflow(log) {
        // set log
        this.log = log

        logHeader()
    }

    // display the title
    def logHeader() {
        log.info workflowHeader
        log.info workflowName
    }
}
