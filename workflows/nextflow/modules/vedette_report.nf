process VEDETTE_REPORT {
    tag "cohort"
    publishDir "${params.outdir}", mode: params.publish_mode, overwrite: true

    input:
    path cohort_vcf
    path cohort_index
    path mendel_report
    path mendel_summary
    path trio_sweep
    path mosaic_tsvs, stageAs: 'mosaic/*'

    output:
    path "report.html", emit: report
    path "families", emit: families, optional: true

    script:
    def ped_arg = params.pedigree ? "--pedigree '${params.pedigree}'" : ""
    def genes_arg = params.genes_gff ? "--genes-gff '${params.genes_gff}'" : ""
    def ndd_arg = params.ndd_tsv ? "--ndd-tsv '${params.ndd_tsv}'" : ""
    def segdup_arg = params.segdup_bed ? "--segdup-bed '${params.segdup_bed}'" : ""
    def lcr_arg = params.lcr_bed ? "--lcr-bed '${params.lcr_bed}'" : ""
    // CLI `--family_workspaces false` reaches here as the STRING "false",
    // which is truthy in Groovy — compare the textual value, not truthiness
    def fam_arg = params.family_workspaces.toString().toLowerCase() in ["false", "0", "no"] \
        ? "--no-family-workspaces" : ""
    def events_arg = params.family_events.toString().toLowerCase() in ["true", "1", "yes"] \
        ? "--family-events" : ""
    def events_pass_arg = params.events_pass_only.toString().toLowerCase() in ["true", "1", "yes"] \
        ? "--events-pass-only" : ""
    def mosaic_arg = params.mosaic.toString().toLowerCase() in ["true", "1", "yes"] \
        ? "--mosaic-dir mosaic" : ""
    def lq_list = params.low_quality_prefixes instanceof List
        ? params.low_quality_prefixes
        : params.low_quality_prefixes.toString().tokenize(',')
    def lq_arg = lq_list.findAll { it?.trim() }.collect { "--low-quality-prefix ${it.trim()}" }.join(' ')
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
    ${params.vedette_bin} report rundir --jobs ${task.cpus} ${ped_arg} ${genes_arg} ${ndd_arg} ${fam_arg} ${events_arg} ${events_pass_arg} ${segdup_arg} ${lcr_arg} ${mosaic_arg} ${lq_arg} \\
        --cohort-name "${params.cohort_name ?: 'cohorte'}" \\
        --out report.html
    if [ -d rundir/families ]; then mv rundir/families families; fi
    """
}
