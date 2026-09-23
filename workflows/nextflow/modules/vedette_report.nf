process VEDETTE_REPORT {
    tag "cohort"
    publishDir "${params.outdir}", mode: params.publish_mode, overwrite: true

    input:
    path cohort_vcf
    path cohort_index
    path mendel_report
    path mendel_summary
    path trio_sweep

    output:
    path "report.html", emit: report
    path "families", emit: families, optional: true

    script:
    def ped_arg = params.pedigree ? "--pedigree '${params.pedigree}'" : ""
    def genes_arg = params.genes_gff ? "--genes-gff '${params.genes_gff}'" : ""
    def ndd_arg = params.ndd_tsv ? "--ndd-tsv '${params.ndd_tsv}'" : ""
    def fam_arg = params.family_workspaces ? "" : "--no-family-workspaces"
    """
    mkdir -p rundir/cohort
    # published per-sample evidence stores: wave-QC column of the sample table
    if [ -d "${params.outdir}/samples" ]; then
        ln -s "${params.outdir}/samples" rundir/samples
    fi
    ln -s \$(readlink -f ${cohort_vcf}) rundir/cohort/cohort.vcf.gz
    ln -s \$(readlink -f ${cohort_index}) rundir/cohort/cohort.vcf.gz.tbi
    ln -s \$(readlink -f ${mendel_report}) rundir/cohort/mendel_report.tsv
    ln -s \$(readlink -f ${mendel_summary}) rundir/cohort/mendel_summary.tsv
    ln -s \$(readlink -f ${trio_sweep}) rundir/cohort/trio_sweep.tsv
    ${params.vedette_bin} report rundir ${ped_arg} ${genes_arg} ${ndd_arg} ${fam_arg} \\
        --cohort-name "${params.cohort_name ?: 'cohorte'}" \\
        --out report.html
    if [ -d rundir/families ]; then mv rundir/families families; fi
    """
}
