/*
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-27
Purpose: Nextflow module
*/

process SamStatsMultiQC {
    tag "${runName}"

    container 'ewels/multiqc:v1.11'

    publishDir "${params.baseDirReport}/align", mode: 'copy', pattern: '*.html'
    publishDir "${params.baseDirData}/align", mode: 'copy', pattern: '*multiqc_data*'

    input:
        file sST
        file sIX
        val runName

    output:
        path '*'

    script:
        """
        multiqc -n ${runName}__fsp_bt2_sST_sIX ${sST} ${sIX}
        """
}
