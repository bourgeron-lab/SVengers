process VEDETTE_CALL {
    tag "$sample"
    publishDir path: { "${params.outdir}/samples/${sample}" }, mode: params.publish_mode, overwrite: true

    input:
    tuple val(sample), val(cram), val(depth_bg), val(vaf_bg)

    output:
    tuple val(sample),
          path("${sample}.vedette.vcf.gz"),
          path("${sample}.vedette.vcf.gz.tbi"),
          path("${sample}.evidence.zarr"), emit: calls

    script:
    def cram_arg = cram ? "--cram '${cram}'" : ""
    def dbg_arg = depth_bg ? "--depth-bedgraph '${depth_bg}'" : ""
    def vbg_arg = vaf_bg ? "--vaf-bedgraph '${vaf_bg}'" : ""
    def pon_arg = params.pon ? "--pon '${params.pon}'" : ""
    def contig_list = params.contigs instanceof List
        ? params.contigs
        : params.contigs.toString().tokenize(',')
    def contig_args = params.contigs ? contig_list.collect { "--contig ${it.trim()}" }.join(' ') : ""
    """
    ${params.vedette_bin} call \\
        --sample ${sample} \\
        --out-dir . \\
        --config ${params.vedette_config} \\
        --resources ${params.resources_dir} \\
        ${cram_arg} ${dbg_arg} ${vbg_arg} ${pon_arg} ${contig_args}
    """
}
