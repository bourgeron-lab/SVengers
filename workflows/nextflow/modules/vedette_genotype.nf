process VEDETTE_GENOTYPE {
    tag "$sample"
    publishDir "${params.outdir}/genotypes", mode: params.publish_mode, overwrite: true

    input:
    tuple val(sample), path(store)
    path sites

    output:
    tuple val(sample), path("${sample}.gt.tsv"), emit: table

    script:
    """
    ${params.vedette_bin} genotype ${sites} \\
        --store ${store} \\
        --out ${sample}.gt.tsv \\
        --config ${params.vedette_config}
    """
}
