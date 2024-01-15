#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { MetaAnalyseCohortsPerGene } from './modules/MetaAnalyseCohorts'
include { PerCohortAnalysisPerGene } from './modules/PerCohortAnalysis'
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
params.variants_percohort = 'NO_FILE'
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
    .splitCsv(header: true, sep: '\t', strip: true).view()

cohort_ch = input_ch.map{row -> row.cohort}.collect()
encoded_ch = input_ch.map{row -> row.encoded}.collect()
genotype_ch = input_ch.map{row -> row.genotype}.collect()
expression_ch = input_ch.map{row -> row.expression}.collect()
partial_derivatives_ch = input_ch.map{row -> row.partial_derivatives}.collect()
snp_inclusion_ch = input_ch.map{row -> row.snp_inclusion}.collect()
gene_inclusion_ch = input_ch.map{row -> row.gene_inclusion}.collect()

all_genes_ch = input_ch
  .map{row -> file(row.gene_inclusion)}
  .splitCsv(header: true)
  .map{gene_row -> gene_row.ID}
  .unique()

genes_per_cohort_ch = Channel.fromPath(params.genes_percohort)
  .splitCsv( header:true ).map { row -> row.ID }

gene_chunk_ch = Channel.of('ID').concat(genes_per_cohort_ch)
  .collectFile(name: 'gene_chunk_ch.txt', keepHeader:false, newLine:true, sort: false, skip:0).view().splitText( by:params.gene_chunk_size, keepHeader:true, file:true ).view()

variants_percohort_ch = Channel.fromPath(params.variants_percohort)
  .collect().view()

mapper = file(params.mapperpath)

// Optional arguments
th = params.th

if (params.covariates) {
  covariate_file = channel.fromPath(params.covariates).collect().view()
} else {covariate_file = Channel.value("")}

workflow {

  PerCohortAnalysisResult = PerCohortAnalysisPerGene(th, gene_chunk_ch, variants_percohort_ch, mapper, covariate_file, cohort_ch, encoded_ch, genotype_ch, expression_ch, partial_derivatives_ch, snp_inclusion_ch, gene_inclusion_ch)

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
