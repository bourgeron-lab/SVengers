process VEDETTE_FILTER {
    tag "$sample"
    publishDir path: { "${params.outdir}/samples/${sample}" }, mode: params.publish_mode, overwrite: true, failOnError: false

    input:
    tuple val(sample), path(vcf), path(tbi), path(store)

    output:
    tuple val(sample),
          path("${sample}.filtered.vcf.gz"),
          path("${sample}.filtered.vcf.gz.tbi"),
          path(store), emit: filtered
    path "${sample}.features.tsv", emit: features

    script:
    def model_arg = params.model_dir ? "--model-dir '${params.model_dir}'" : ""
    def artifact_arg = params.artifact_bed ? "--artifact-bed '${params.artifact_bed}'" : ""
    """
    ${params.vedette_bin} filter ${vcf} \\
        --store ${store} \\
        --out ${sample}.filtered.vcf.gz \\
        --config ${params.vedette_config} \\
        --features-out ${sample}.features.tsv \\
        ${model_arg} ${artifact_arg}
    """
}
