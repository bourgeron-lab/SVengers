process SVENGERS_BND {
    tag "cohort"
    // per-sample outputs of a run-level step: published beside each sample's calls
    publishDir path: "${params.outdir}", mode: params.publish_mode, overwrite: true, failOnError: false,
               saveAs: { fn -> fn.startsWith('samples/') ? fn : null }

    input:
    path filtered, stageAs: 'vcf/*'
    path stores, stageAs: 'store/*'
    path mosaic_tsvs, stageAs: 'mosaic/*'

    output:
    path "samples/*/*.{bnd.tsv,bnd.vcf.gz,bnd.vcf.gz.tbi,terminal.tsv}", emit: files

    script:
    """
    # a run directory of links: samples/<S>/<S>.{evidence.zarr,filtered.vcf.gz,mosaic.tsv}
    for st in store/*.evidence.zarr; do
        s=\$(basename \$st .evidence.zarr)
        mkdir -p samples/\$s
        ln -s \$(readlink -f \$st) samples/\$s/\$s.evidence.zarr
        for f in vcf/\$s.filtered.vcf.gz vcf/\$s.filtered.vcf.gz.tbi mosaic/\$s.mosaic.tsv; do
            if [ -e \$f ]; then ln -s \$(readlink -f \$f) samples/\$s/\$(basename \$f); fi
        done
    done
    ${params.svengers_bin} bnd --run-dir . --fasta ${params.fasta} --jobs ${task.cpus}
    """
}
