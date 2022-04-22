include { DeepToolsMultiBamSummary } from "${projectDir}/modules/align/DeepToolsMultiBamSummary.nf"

workflow DeepToolsMultiBamSWF {
    take:
        bams
        bais
        outFileName
    
    main:
        DeepToolsMultiBamSummary(
            bams,
            bais,
            outFileName
        )
    
    emit:
        corMatrix = DeepToolsMultiBamSummary.out.corMatrix
        PCAMatrix = DeepToolsMultiBamSummary.out.PCAMatrix
}
