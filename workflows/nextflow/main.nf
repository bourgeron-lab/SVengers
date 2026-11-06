#!/usr/bin/env nextflow

/*
========================================================================================
    Vedette cohort CNV-calling pipeline
========================================================================================
    call -> filter -> merge -> genotype -> cohort -> ped -> (mosaic, bnd) -> report

    Input: a Vedette manifest TSV (sample_id, cram, batch, sex, depth_bedgraph,
    vaf_bedgraph) built with `vedette manifest discover`, plus optionally a
    pedigree. Sample-level steps run in parallel; cohort steps fan in.
    Restartable with -resume.
----------------------------------------------------------------------------------------
*/

include { VEDETTE_CALL } from './modules/vedette_call'
include { VEDETTE_FILTER } from './modules/vedette_filter'
include { VEDETTE_MERGE } from './modules/vedette_merge'
include { VEDETTE_GENOTYPE } from './modules/vedette_genotype'
include { VEDETTE_COHORT_FILTER } from './modules/vedette_cohort_filter'
include { VEDETTE_PED } from './modules/vedette_ped'
include { VEDETTE_MOSAIC } from './modules/vedette_mosaic'
include { VEDETTE_BND } from './modules/vedette_bnd'
include { VEDETTE_REPORT } from './modules/vedette_report'

def validatedSteps() {
    def valid_steps = ['call', 'filter', 'merge', 'genotype', 'cohort', 'ped', 'report']
    if (!params.manifest) {
        error "ERROR: --manifest parameter is required (vedette manifest discover ...)"
    }
    // --steps on the CLI arrives as a comma-separated string; params files give a list
    def steps = (params.steps instanceof List ? params.steps : params.steps.toString().tokenize(','))
        .collect { it.trim() }
    if (!steps) {
        error "ERROR: --steps parameter is required. Available steps: ${valid_steps.join(', ')}"
    }
    def invalid = steps - valid_steps
    if (invalid) {
        error "ERROR: Invalid steps: ${invalid.join(', ')}. Valid steps: ${valid_steps.join(', ')}"
    }
    if ('ped' in steps && !params.pedigree) {
        error "ERROR: step 'ped' requires --pedigree"
    }
    // downstream steps need their upstream in the same run, with one exception:
    // starting at 'filter' discovers published calls (--calls_dir), the re-filter
    // use case after a model update
    def wanted = valid_steps.findAll { it in steps }
    if (wanted != valid_steps[0..wanted.size() - 1]
            && !(wanted == valid_steps[1..wanted.size()])) {
        error "ERROR: steps must be a prefix of ${valid_steps.join(' -> ')} (or start at 'filter' with --calls_dir; got: ${wanted.join(', ')})"
    }
    if (wanted[0] == 'filter' && !params.calls_dir) {
        error "ERROR: starting at 'filter' requires --calls_dir (published <sample>/<sample>.vedette.vcf.gz + evidence store)"
    }
    return wanted
}

// Re-filter mode: the report (and a later `vedette report` by hand) reads each
// store at <outdir>/samples/<S>/, but only the filter outputs are published
// there -- a 2670-sample re-filter lost the sex-chromosome section and wave
// alerts without a word (2026-11-02). Symlink the published call instead: no
// copy of the stores, and no publishDir on a passthrough input (the symlink
// cycle of VEDETTE_FILTER). Entries already present are left alone.
def linkPublishedCall(String sample, src_dir) {
    def dest_dir = file("${params.outdir}/samples/${sample}")
    if (dest_dir.exists() && java.nio.file.Files.isSameFile(dest_dir, src_dir)) {
        return  // --calls_dir is this run's own samples/
    }
    dest_dir.mkdirs()
    ["${sample}.evidence.zarr", "${sample}.vedette.vcf.gz", "${sample}.vedette.vcf.gz.tbi"].each { fn ->
        def src = src_dir.resolve(fn)
        def dest = dest_dir.resolve(fn)
        if (src.exists() && !java.nio.file.Files.exists(dest, java.nio.file.LinkOption.NOFOLLOW_LINKS)) {
            java.nio.file.Files.createSymbolicLink(dest, src.toRealPath())
        }
    }
}

workflow {

    def steps = validatedSteps()

    log.info """
    ========================================================================================
                              Vedette CNV cohort pipeline
    ========================================================================================
    Manifest         : ${params.manifest}
    Pedigree         : ${params.pedigree ?: '(none)'}
    Steps            : ${steps.join(', ')}
    Vedette config   : ${params.vedette_config}
    PoN              : ${params.pon ?: '(none)'}
    Filter model     : ${params.model_dir ?: '(tier-1 only)'}
    Output directory : ${params.outdir}
    Work directory   : ${workflow.workDir}
    ========================================================================================
    """

    // manifest rows: sample_id, cram, batch, sex, depth_bedgraph, vaf_bedgraph
    def samples = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            [
                row.sample_id,
                row.cram,
                row.depth_bedgraph ?: '',
                row.vaf_bedgraph ?: '',
            ]
        }

    def calls = Channel.empty()
    if ('call' in steps) {
        VEDETTE_CALL(samples)
        calls = VEDETTE_CALL.out.calls
    } else {
        // re-filter mode: pick up published calls for the manifest's samples
        calls = samples.map { sample, cram, dbg, vbg ->
            def d = file("${params.calls_dir}/${sample}")
            def vcf = file("${d}/${sample}.vedette.vcf.gz")
            def store = file("${d}/${sample}.evidence.zarr")
            if (!vcf.exists() || !store.exists()) {
                error "ERROR: published call missing for ${sample} under ${params.calls_dir}"
            }
            linkPublishedCall(sample, d)
            tuple(sample, vcf, file("${vcf}.tbi"), store)
        }
    }

    def filtered = Channel.empty()
    def sites = Channel.empty()
    def genotype_tables = Channel.empty()

    if ('filter' in steps) {
        VEDETTE_FILTER(calls)
        filtered = VEDETTE_FILTER.out.filtered
    }

    if ('merge' in steps) {
        VEDETTE_MERGE(filtered.map { sample, vcf, tbi, store -> vcf }.collect())
        sites = VEDETTE_MERGE.out.sites
    }

    if ('genotype' in steps) {
        VEDETTE_GENOTYPE(
            filtered.map { sample, vcf, tbi, store -> [sample, store] },
            sites,
        )
        genotype_tables = VEDETTE_GENOTYPE.out.table
    }

    if ('cohort' in steps) {
        VEDETTE_COHORT_FILTER(
            sites,
            genotype_tables.map { sample, tsv -> tsv }.collect(),
        )
    }

    if ('ped' in steps) {
        VEDETTE_PED(VEDETTE_COHORT_FILTER.out.cohort_vcf, VEDETTE_COHORT_FILTER.out.cohort_index)
    }

    // chromosomal-scale mosaics, per sample from its store; the report stages
    // the TSVs (an input, not the published copies, so it never races them)
    def mosaic_tsvs = Channel.value([])
    if ('report' in steps && params.mosaic.toString().toLowerCase() in ["true", "1", "yes"]) {
        VEDETTE_MOSAIC(calls.map { sample, vcf, tbi, store -> [sample, store] })
        mosaic_tsvs = VEDETTE_MOSAIC.out.tsv.collect().ifEmpty([])
    }

    // translocation / inversion junctions and the mechanism of terminal losses,
    // over the whole run (cohort flags); needs the filtered calls of this run
    def bnd_files = Channel.value([])
    if ('report' in steps && 'filter' in steps
            && params.bnd.toString().toLowerCase() in ["true", "1", "yes"]) {
        VEDETTE_BND(
            filtered.map { sample, vcf, tbi, store -> [vcf, tbi] }.flatten().collect(),
            filtered.map { sample, vcf, tbi, store -> store }.collect(),
            mosaic_tsvs,
        )
        bnd_files = VEDETTE_BND.out.files.collect().ifEmpty([])
    }

    if ('report' in steps) {
        VEDETTE_REPORT(
            VEDETTE_COHORT_FILTER.out.cohort_vcf,
            VEDETTE_COHORT_FILTER.out.cohort_index,
            VEDETTE_PED.out.report,
            VEDETTE_PED.out.summary,
            VEDETTE_PED.out.sweep,
            mosaic_tsvs,
            bnd_files,
        )
    }
}
