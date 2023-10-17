/**
 * 'Quiet grep' a file path for a query value.
 * 
 * Check if a file contains a line with the query.
 * This should essentially mimic running `grep -q` at the command line.
 * The goal here is speed, so as soon as a query is found the iteration is broken and a `true` value is returned.
 *
 * @params path A file path.
 * @params query A query value.
 *
 * @return boolean Is the query is found in path.
 */
static def grepq(path, query) {
    def containsQuery = false
    
    // iterate through path line by line.
    // set that path contains query and exit the iteration as soon as query is found
    path.eachLine{ line ->
        if(line.contains(query)) {
            containsQuery = true
            return // exit the eachLine closure
        }
    }

    return containsQuery
}
