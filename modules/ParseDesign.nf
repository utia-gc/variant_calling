process ParseDesign {
    tag "${design}"

    container "kyclark/tiny_python_projects:0.2.0"

    publishDir "${params.baseDirData}/design", mode: 'copy'

    input:
        path design

    output:
        path '*.csv', emit: csv

    script:
    """
    parse_design.py \
        ${design}
    """
}
