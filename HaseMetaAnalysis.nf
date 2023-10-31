#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { MetaAnalyseCohortsPerGene } from './modules/MetaAnalyseCohorts'
include { PerCohortAnalysis } from './modules/PerCohortAnalysis'
include { SubsetGenesInclusion } from './modules/SubsetGenesInclusion'
include { Partition; Partition as PartitionPerCohort; CleanPartition } from './modules/Partition'
include { ListPhenotypes } from './modules/ListPhenotypes'
include { Combine; Combine as CombinePerCohort; CleanCombine } from './modules/Combine'

def helpmessage() {

log.info"""

HASE meta-analyzer ~ v${workflow.manifest.version}"
=================================================
Pipeline for running encoded meta-analysis over multiple studies in the central site.

Pay attention that it does not explicitly overwrite the output folder, so clean it before re-running the pipeline.

Usage:

nextflow run HaseMetaAnalysis.nf \
--mastertable [path to master table file] \
--mapperpath [path to mapper files] \
--outdir [path to results folder]

Mandatory arguments:
--mastertable     Path to the master file with following columns. cohort: cohort name, genotype: path to encoded genotype data, expression: path to encoded expression data, partial_derivatives: path to partial derivatives folder, encoded: encoded status 0 or 1, snp_inclusion: path to SNPs inclusion file (header "ID", and SNP IDs on each row) and gene_inclusion: path to gene inclusion file (header "ID", and gene IDs on each row).
--mapperpath      Path to the mappers folder, containing the mapper files for each cohort.
--outdir      Path to results folder where meta-analysis outputs are written.
--covariates      What covariates to include into the analysis for each cohort. If specified, the it has to be the file which has each line formatted as following: [cohort name] [index of covariate to include 1] [index of covariate to include 2] etc. Indexing is 1-based and defaults to all covariates.


Optional arguments:
--th             T-statistic threshold for writing out the results. Defaults to 0.
--genes_percohort File with genes to output cohort specific results for.
--chunks         In how many chunks to run the analysis. Defaults to 100.

""".stripIndent()

}

//Default parameters
params.mastertable = ''
params.genes_percohort = ''
params.gene_filter = ''
params.outdir = ''
params.mapperpath = ''
params.th = 0
params.covariates = ''
params.gene_chunk_size = 100

log.info """=================================================
HASE meta-analyzer v${workflow.manifest.version}"
================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Master table']                             = params.mastertable
summary['Genes per cohort']                         = params.genes_percohort
summary['Input mapper path']                        = params.mapperpath
summary['Output directory']                         = params.outdir
summary['T-statistic threshold']                    = params.th
summary['Covariate filter']                         = params.covariates
summary['Nr. of chunks']                            = params.chunks

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "================================================="


// Process input file paths

input_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)

cohort_ch = input_ch.map{row -> row.cohort}.collect()
encoded_ch = input_ch.map{row -> row.encoded}.collect()
genotype_ch = input_ch.map{row -> row.genotype}.collect()
expression_ch = input_ch.map{row -> row.expression}.collect()
partial_derivatives_ch = input_ch.map{row -> row.partial_derivatives}.collect()
snp_inclusion_ch = input_ch.map{row -> row.snp_inclusion}.collect()
gene_inclusion_ch = input_ch.map{row -> row.gene_inclusion}.collect()

if (params.gene_filter) {
  gene_chunk_ch = Channel.fromPath(params.gene_filter)
   .ifEmpty { error "Cannot find gene filter from: ${params.gene_filter}" }
   .splitText( by:params.gene_chunk_size, keepHeader:true, file:true)

} else {
  all_genes_ch = input_ch.map{row -> row.gene_inclusion}
   .splitCsv(header: true)
   .map{row -> row.ID}
   .flatten()
   .unique()

  gene_chunk_ch = Channel.of('ID').concat(all_genes_ch)
  .collectFile()
  .splitText( by:params.gene_chunk_size, keepHeader:true, file:true)
}

if (params.genes_percohort) {
  gene_percohort_ch = Channel.fromPath(params.genes_percohort)
   .ifEmpty { error "Cannot find gene per cohort from: ${params.genes_percohort}" }
}

if (params.variants_percohort) {
  variants_percohort_ch = Channel.fromPath(params.variants_percohort)
   .ifEmpty { error "Cannot find variant per cohort from: ${params.variants_percohort}" }
}

mapper = file(params.mapperpath)

test_ch = Channel.fromPath(params.mastertable)

// variant_reference_ch = Channel.fromPath(params.reference).collect()

// Optional arguments
th = params.th
chunks = params.chunks

if (params.covariates) {
  covariate_file = file(params.covariates)
} else {covariate_file = Channel.value("")}

workflow {

  if (params.genes_percohort) {
    PerCohortAnalysisResult = PerCohortAnalysis(th, gene_chunk_ch, gene_percohort_ch, mapper, covariate_file, cohort_ch, encoded_ch, genotype_ch, expression_ch, partial_derivatives_ch, snp_inclusion_ch, gene_inclusion_ch)
  }

  else if (params.gene_filter) {
    //MetaAnalysisResult = MetaAnalyseCohortsPerGene(th, gene_chunk_ch, mapper, covariate_file, input_ch)

  }

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
