class Args {
    LinkedHashMap<String, String> defaults = new LinkedHashMap<String, String>()
    LinkedHashMap<String, String> dynamics = new LinkedHashMap<String, String>()
    LinkedHashMap<String, String> users    = new LinkedHashMap<String, String>()

    public String buildArgsString() {
        LinkedHashMap consensusArgs = mergeArgs()

        def builder = new StringBuilder()
        consensusArgs.each { parameter, argument ->
            // separate key:value pairs by spaces
            if (builder.length() > 0) {
                builder.append(' ')
            }
            argument ? builder.append("${parameter} ${argument}") : builder.append("${parameter}")
        }

        return builder.toString()
    }

    private LinkedHashMap mergeArgs() {
        LinkedHashMap mergedArgs = [:]
        mergedArgs += defaults
        mergedArgs += dynamics
        mergedArgs += users

        return mergedArgs
    }
}

// public static String buildArgsString(LinkedHashMap args) {
//     def builder = new StringBuilder()
//     args.each { parameter, argument ->
//         // separate key:value pairs by spaces
//         if (builder.length() > 0) {
//             builder.append(' ')
//         }
//         argument ? builder.append("${parameter} ${argument}") : builder.append("${parameter}")
//     }

//     return builder.toString()
// }
// public static LinkedHashMap deepMerge(LinkedHashMap... maps) {
//     def mergedMap = [:]
//     maps.each { map ->
//         mergedMap += [ *: map ]
//     }

//     return mergedMap
// }

// LinkedHashMap<String, String> argsDefault = [
//     '--param1': 'arg1',
//     '--flag1': ''
// ]
// LinkedHashMap<String, String> argsDynamic = [
//     '--paired': ''
// ]
// LinkedHashMap<String, String> argsUser = [
//     '--param1': 'different_arg',
//     '--flag1': ''
// ]

// def args = new Args()

// println(args)
