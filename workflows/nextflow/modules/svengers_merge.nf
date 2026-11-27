process SVENGERS_MERGE {
    tag "cohort"
    publishDir "${params.outdir}/cohort", mode: params.publish_mode, overwrite: true

    input:
    path filtered_vcfs

    output:
    path "sites.vcf.gz", emit: sites
    path "sites.vcf.gz.tbi", emit: sites_index
    path "sites.members.tsv.gz", emit: members

    script:
    """
    ls -1 ${filtered_vcfs} > vcf_list.txt
    ${params.svengers_bin} merge @vcf_list.txt \\
        --jobs ${task.cpus} \\
        --out sites.vcf.gz \\
        --fasta ${params.fasta}
    """
}
