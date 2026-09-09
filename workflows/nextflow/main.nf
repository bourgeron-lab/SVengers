#!/usr/bin/env nextflow

/*
========================================================================================
    Vedette cohort CNV-calling pipeline
========================================================================================
    call -> filter -> merge -> genotype -> cohort -> ped

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
    // downstream steps need their upstream in the same run (no cross-run discovery in v1)
    def wanted = valid_steps.findAll { it in steps }
    if (wanted != valid_steps[0..wanted.size() - 1]) {
        error "ERROR: steps must be a prefix of ${valid_steps.join(' -> ')} (got: ${wanted.join(', ')})"
    }
    return wanted
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

    VEDETTE_CALL(samples)

    def filtered = Channel.empty()
    def sites = Channel.empty()
    def genotype_tables = Channel.empty()

    if ('filter' in steps) {
        VEDETTE_FILTER(VEDETTE_CALL.out.calls)
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
        VEDETTE_PED(VEDETTE_COHORT_FILTER.out.cohort_vcf)
    }

    if ('report' in steps) {
        VEDETTE_REPORT(
            VEDETTE_COHORT_FILTER.out.cohort_vcf,
            VEDETTE_COHORT_FILTER.out.cohort_index,
            VEDETTE_PED.out.report,
            VEDETTE_PED.out.summary,
            VEDETTE_PED.out.sweep,
        )
    }
}
