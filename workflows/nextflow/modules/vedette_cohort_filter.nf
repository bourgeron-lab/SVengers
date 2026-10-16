process VEDETTE_COHORT_FILTER {
    tag "cohort"
    publishDir "${params.outdir}/cohort", mode: params.publish_mode, overwrite: true

    input:
    path sites
    path genotype_tables

    output:
    path "cohort.vcf.gz", emit: cohort_vcf
    path "cohort.vcf.gz.tbi", emit: cohort_index

    script:
    def ped_arg = params.pedigree ? "--pedigree '${params.pedigree}'" : ""
    def known_arg = params.known_sites ? "--known-sites '${params.known_sites}'" : ""
    """
    ls -1 ${genotype_tables} > gt_list.txt
    ${params.vedette_bin} cohort-filter ${sites} \\
        --jobs ${task.cpus} \\
        --genotypes @gt_list.txt \\
        --out cohort.vcf.gz \\
        ${ped_arg} ${known_arg}
    """
}
