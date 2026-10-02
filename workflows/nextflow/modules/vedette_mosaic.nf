process VEDETTE_MOSAIC {
    tag "$sample"
    publishDir path: { "${params.outdir}/samples/${sample}" }, mode: params.publish_mode, overwrite: true, failOnError: false

    input:
    tuple val(sample), path(store)

    output:
    path "${sample}.mosaic.tsv", emit: tsv

    script:
    """
    ${params.vedette_bin} mosaic ${store} --out ${sample}.mosaic.tsv
    """
}
