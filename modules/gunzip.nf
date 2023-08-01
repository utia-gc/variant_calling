process GUNZIP {
    tag "${archive}"

    label 'base'

    label 'lil_mem'

    input:
    path archive
    
    output:
    path "${archive.toString() - '.gz'}", emit: gunzip
    
    script:
    """
    gunzip \
        --force \
        ${archive}
    """
}
