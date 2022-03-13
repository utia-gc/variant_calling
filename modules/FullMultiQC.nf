/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-03-13
Purpose: multiqc report for full pipeline
*/

process FullMultiQC {
    tag "${runName}"

    container 'ewels/multiqc:v1.11'

    publishDir "${params.baseDirReport}", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}",   mode: 'copy', pattern: '*multiqc_data*'

    input:
        val  runName
        file config
        file rawFastqc
        file trimFastqc
        file alignGenomeStats
        file alignGenomeIdxstats

    output:
        path "*"

    script:
        """
        multiqc \
            -n ${runName} -i ${runName} \
            -c ${config} \
            ${rawFastqc} \
            ${trimFastqc} \
            ${alignGenomeStats} \
            ${alignGenomeIdxstats}
        """
}
