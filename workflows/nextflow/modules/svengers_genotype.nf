process SVENGERS_GENOTYPE {
    tag "$sample"
    publishDir "${params.outdir}/genotypes", mode: params.publish_mode, overwrite: true

    input:
    tuple val(sample), path(store)
    path sites
    path members

    output:
    tuple val(sample), path("${sample}.gt.tsv"), emit: table

    script:
    """
    ${params.svengers_bin} genotype ${sites} \\
        --store ${store} \\
        --out ${sample}.gt.tsv \\
        --members ${members} \\
        --config ${params.svengers_config}
    """
}
