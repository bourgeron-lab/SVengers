#!/usr/bin/env nextflow

/*
========================================================================================
    SVengers cohort CNV-calling pipeline
========================================================================================
    call -> filter -> merge -> genotype -> cohort -> ped -> (mosaic, bnd) -> report

    Input: an SVengers manifest TSV (sample_id, cram, batch, sex, depth_bedgraph,
    vaf_bedgraph) built with `svengers manifest discover`, plus optionally a
    pedigree. Sample-level steps run in parallel; cohort steps fan in.
    Restartable with -resume.
----------------------------------------------------------------------------------------
*/

include { SVENGERS_CALL } from './modules/svengers_call'
include { SVENGERS_FILTER } from './modules/svengers_filter'
include { SVENGERS_MERGE } from './modules/svengers_merge'
include { SVENGERS_GENOTYPE } from './modules/svengers_genotype'
include { SVENGERS_COHORT_FILTER } from './modules/svengers_cohort_filter'
include { SVENGERS_PED } from './modules/svengers_ped'
include { SVENGERS_MOSAIC } from './modules/svengers_mosaic'
include { SVENGERS_BND } from './modules/svengers_bnd'
include { SVENGERS_REPORT } from './modules/svengers_report'

def validatedSteps() {
    def valid_steps = ['call', 'filter', 'merge', 'genotype', 'cohort', 'ped', 'report']
    if (!params.manifest) {
        error "ERROR: --manifest parameter is required (svengers manifest discover ...)"
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
        error "ERROR: starting at 'filter' requires --calls_dir (published <sample>/<sample>.svengers.vcf.gz, or .vedette.vcf.gz before v2.0, + evidence store)"
    }
    return wanted
}

// A sample's published raw calls: <S>.svengers.vcf.gz, or <S>.vedette.vcf.gz
// for calls made before v2.0 (when SVengers was Vedette)
def publishedRawCalls(String sample, dir) {
    def vcf = file("${dir}/${sample}.svengers.vcf.gz")
    def legacy = file("${dir}/${sample}.vedette.vcf.gz")
    return (!vcf.exists() && legacy.exists()) ? legacy : vcf
}

// Re-filter mode: the report (and a later `svengers report` by hand) reads each
// store at <outdir>/samples/<S>/, but only the filter outputs are published
// there -- a 2670-sample re-filter lost the sex-chromosome section and wave
// alerts without a word (2026-11-02). Symlink the published call instead: no
// copy of the stores, and no publishDir on a passthrough input (the symlink
// cycle of SVENGERS_FILTER). Entries already present are left alone.
def linkPublishedCall(String sample, src_dir) {
    def dest_dir = file("${params.outdir}/samples/${sample}")
    if (dest_dir.exists() && java.nio.file.Files.isSameFile(dest_dir, src_dir)) {
        return  // --calls_dir is this run's own samples/
    }
    dest_dir.mkdirs()
    def vcf_name = publishedRawCalls(sample, src_dir).getName()
    ["${sample}.evidence.zarr", vcf_name, "${vcf_name}.tbi"].each { fn ->
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
                              SVengers CNV cohort pipeline
    ========================================================================================
    Manifest         : ${params.manifest}
    Pedigree         : ${params.pedigree ?: '(none)'}
    Steps            : ${steps.join(', ')}
    SVengers config  : ${params.svengers_config}
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
        SVENGERS_CALL(samples)
        calls = SVENGERS_CALL.out.calls
    } else {
        // re-filter mode: pick up published calls for the manifest's samples
        calls = samples.map { sample, cram, dbg, vbg ->
            def d = file("${params.calls_dir}/${sample}")
            def vcf = publishedRawCalls(sample, d)
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
        SVENGERS_FILTER(calls)
        filtered = SVENGERS_FILTER.out.filtered
    }

    if ('merge' in steps) {
        SVENGERS_MERGE(filtered.map { sample, vcf, tbi, store -> vcf }.collect())
        sites = SVENGERS_MERGE.out.sites
    }

    if ('genotype' in steps) {
        SVENGERS_GENOTYPE(
            filtered.map { sample, vcf, tbi, store -> [sample, store] },
            sites,
        )
        genotype_tables = SVENGERS_GENOTYPE.out.table
    }

    if ('cohort' in steps) {
        SVENGERS_COHORT_FILTER(
            sites,
            genotype_tables.map { sample, tsv -> tsv }.collect(),
        )
    }

    if ('ped' in steps) {
        SVENGERS_PED(SVENGERS_COHORT_FILTER.out.cohort_vcf, SVENGERS_COHORT_FILTER.out.cohort_index)
    }

    // chromosomal-scale mosaics, per sample from its store; the report stages
    // the TSVs (an input, not the published copies, so it never races them)
    def mosaic_tsvs = Channel.value([])
    if ('report' in steps && params.mosaic.toString().toLowerCase() in ["true", "1", "yes"]) {
        SVENGERS_MOSAIC(calls.map { sample, vcf, tbi, store -> [sample, store] })
        mosaic_tsvs = SVENGERS_MOSAIC.out.tsv.collect().ifEmpty([])
    }

    // translocation / inversion junctions and the mechanism of terminal losses,
    // over the whole run (cohort flags); needs the filtered calls of this run
    def bnd_files = Channel.value([])
    if ('report' in steps && 'filter' in steps
            && params.bnd.toString().toLowerCase() in ["true", "1", "yes"]) {
        SVENGERS_BND(
            filtered.map { sample, vcf, tbi, store -> [vcf, tbi] }.flatten().collect(),
            filtered.map { sample, vcf, tbi, store -> store }.collect(),
            mosaic_tsvs,
        )
        bnd_files = SVENGERS_BND.out.files.collect().ifEmpty([])
    }

    if ('report' in steps) {
        SVENGERS_REPORT(
            SVENGERS_COHORT_FILTER.out.cohort_vcf,
            SVENGERS_COHORT_FILTER.out.cohort_index,
            SVENGERS_PED.out.report,
            SVENGERS_PED.out.summary,
            SVENGERS_PED.out.sweep,
            mosaic_tsvs,
            bnd_files,
        )
    }
}
