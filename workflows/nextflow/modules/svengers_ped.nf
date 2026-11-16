process SVENGERS_PED {
    tag "cohort"
    publishDir "${params.outdir}/cohort", mode: params.publish_mode, overwrite: true

    input:
    path cohort_vcf
    // staged next to the VCF: without it the region split falls back to one
    // core (1 h 08 instead of ~9 min on 2670 samples)
    path cohort_index

    output:
    path "mendel_report.tsv", emit: report
    path "mendel_summary.tsv", emit: summary
    path "trio_sweep.tsv", emit: sweep
    path "trio_sweep.png", emit: sweep_plot, optional: true

    script:
    // ClinGen dosage map: a de novo at a dosage-sensitive locus may recur across families
    def clingen_arg = (params.clingen_genes ? "--clingen-genes '${params.clingen_genes}' " : "") \
        + (params.clingen_regions ? "--clingen-regions '${params.clingen_regions}'" : "")
    """
    ${params.svengers_bin} ped-annotate ${cohort_vcf} ${clingen_arg} \\
        --jobs ${task.cpus} \\
        --pedigree '${params.pedigree}' \\
        --out-report mendel_report.tsv \\
        --out-summary mendel_summary.tsv \\
        --calls-dir '${params.outdir}/samples'
    ${params.svengers_bin} bench trios ${cohort_vcf} \\
        --jobs ${task.cpus} \\
        --pedigree '${params.pedigree}' \\
        --out trio_sweep.tsv --plot || true
    # a cohort without complete trios yields no sweep: emit an empty one so
    # downstream report wiring always fires
    touch trio_sweep.tsv
    """
}
