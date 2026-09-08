process VEDETTE_PED {
    tag "cohort"
    publishDir "${params.outdir}/cohort", mode: params.publish_mode, overwrite: true

    input:
    path cohort_vcf

    output:
    path "mendel_report.tsv", emit: report
    path "mendel_summary.tsv", emit: summary
    path "trio_sweep.tsv", emit: sweep, optional: true
    path "trio_sweep.png", emit: sweep_plot, optional: true

    script:
    """
    ${params.vedette_bin} ped-annotate ${cohort_vcf} \\
        --pedigree '${params.pedigree}' \\
        --out-report mendel_report.tsv \\
        --out-summary mendel_summary.tsv
    ${params.vedette_bin} bench trios ${cohort_vcf} \\
        --pedigree '${params.pedigree}' \\
        --out trio_sweep.tsv --plot || true
    """
}
