process FASTQC {
    label 'fastqc'
    label 'lil_mem'

    container = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    publishDir(path: "${publish_dir}/qc/${outdir_name}", mode: "copy")

    input:
        tuple val(metadata), path(reads)
        val outdir_name

    output:
        path("fastqc_${metadata.sampleName}_logs"), emit: fastq_ch

    script:
        if (metadata.readType == 'single') {
            """
            mkdir fastqc_${metadata.sampleName}_logs
            fastqc -o fastqc_${metadata.sampleName}_logs -q ${reads}
            """
        } else {
            """
            mkdir fastqc_${metadata.sampleName}_logs
            fastqc -o fastqc_${metadata.sampleName}_logs -q ${reads}
            """
        }
}
