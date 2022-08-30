#!/usr/bin/env nextflow

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
--outdir          Path to results folder where meta-analysis outputs are written.
--covariates      What covariates to include into the analysis for each cohort. If specified, the it has to be the file which has each line formatted as following: [cohort name] [index of covariate to include 1] [index of covariate to include 2] etc. Indexing is 1-based and defaults to all covariates.


Optional arguments:
--th             T-statistic threshold for writing out the results. Defaults to 0.
--chunks         In how many chunks to run the analysis. Defaults to 100.

""".stripIndent()

}

//Default parameters
params.mastertable = ''
params.outdir = ''
params.mapperpath = ''
params.th = 0
params.covariates = ''
params.chunks = 100
params.permuted = 'no'

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
summary['Input mapper path']                        = params.mapperpath
summary['Output directory']                         = params.outdir
summary['T-statistic threshold']                    = params.th
summary['Covariate filter']                         = params.covariates
summary['Nr. of chunks']                            = params.chunks

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "================================================="

// Count the number of input datasets

// Process input file paths
mastertable = file(params.mastertable)

cohort_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "${row.cohort}" ]}
    .view{}
    
cohort = Channel.value("${cohort_ch}")

genotype_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "${row.genotype}" ]}
    
genotype = Channel.value("${genotype_ch}")

expression_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "${row.expression}" ]}
    .view()
    
expression = Channel.value("${expression_ch}")

partial_derivatives_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "${row.partial_derivatives}" ]}
    
partial_derivatives = Channel.value("${partial_derivatives_ch}")

encoded_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "${row.encoded}" ]}
    
encoded = Channel.value("${encoded_ch}")

snp_inclusion_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "-snp_id_inc ${row.snp_inclusion}" ]}
    
snp_inclusion = Channel.value("${snp_inclusion_ch}")

gene_inclusion_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ "${row.gene_inclusion}" ]}
    
gene_inclusion = Channel.value("-ph_id_inc ${gene_inclusion_ch}")
/* 
mapper = file(params.mapperpath)

// Optional arguments
th = params.th
chunks = params.chunks

if (params.covariates) { 
  covariate_file = file(params.covariates) 
  covariate_filtering = Channel.value("-ci ${covariate_file}")
} else {covariate_filtering = Channel.value("")}

chunk = Channel.from(1..chunks)

process MetaAnalyseCohorts {

    tag "chunk $chunk"

    publishDir "${params.outdir}", mode: 'copy'

    input:
      val chunk from chunk
      path mapper from mapper
      val(th) from th
      val(nr_chunks) from chunks
      val(snp_inclusion) from snp_inclusion
      val(gene_inclusion) from gene_inclusion
      val(covariate_filtering) from covariate_filtering
      val(genotype) from genotype
      val(expression) from expression
      val(partial_derivatives) from partial_derivatives
      val(cohort) from cohort
      val(encoded) from encoded

    output:
      path 'MetaAnalysisResultsEncoded/*.parquet' into MetaAnalysisResultsEncoded

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    python2 !{baseDir}/bin/hase/hase.py \
    -study_name !{cohort} \
    -g !{genotype} \
    -ph !{expression} \
    -derivatives !{partial_derivatives} \
    -mapper !{mapper}/ \
    -o MetaAnalysisResultsEncoded \
    -mode meta-classic \
    -encoded ${encodedstatus} \
    -allow_missingness \
    -node !{nr_chunks} !{chunk} \
    -th !{th} \
    -mapper_chunk 500 \
    -ref_name 1000G-30x_ref \
    -cluster "y" \
    !{snp_inclusion} \
    !{gene_inclusion} \
    !{covariate_filtering}
    '''
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
 */