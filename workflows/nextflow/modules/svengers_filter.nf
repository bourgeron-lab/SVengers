process SVENGERS_FILTER {
    tag "$sample"
    // The evidence store is passed through to the downstream steps, not produced
    // here: publishing it would overwrite the published entry the task read from
    // with a link to the work-dir copy of itself -- a symlink cycle that makes
    // the store unreadable. Publish the filter's own outputs only.
    publishDir path: { "${params.outdir}/samples/${sample}" }, mode: params.publish_mode,
               overwrite: true, failOnError: false,
               saveAs: { fn -> fn.endsWith('.evidence.zarr') ? null : fn }

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
    ${params.svengers_bin} filter ${vcf} \\
        --store ${store} \\
        --out ${sample}.filtered.vcf.gz \\
        --config ${params.svengers_config} \\
        --features-out ${sample}.features.tsv \\
        ${model_arg} ${artifact_arg}
    """
}
