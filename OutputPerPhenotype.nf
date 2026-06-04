#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { OldPartitionPerCohort } from './modules/Partition'
include { ListPhenotypes } from './modules/ListPhenotypes'
include { Combine } from './modules/Combine'

def helpmessage() {

log.info"""

HASE meta-analyzer ~ v${workflow.manifest.version}"
=================================================
Pipeline for running encoded meta-analysis over multiple studies in the central site.

Pay attention that it does not explicitly overwrite the output folder, so clean it before re-running the pipeline.

""".stripIndent()

}

//Default parameters
params.input = ''
params.outdir = ''
params.chunks = 100

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
summary['Input path']                               = params.input
summary['Output directory']                         = params.outdir

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "================================================="


// Process input file paths

parquet = Channel.fromPath(params.input).collect()
variant_reference_ch = Channel.fromPath(params.reference).collect()

chunks = params.chunks
chunk = Channel.from(1..chunks)

workflow {

OldPartitionPerCohort(parquet, chunk)
ListPhenotypes(OldPartitionPerCohort.out.phenotypes.collect())
Combine(OldPartitionPerCohort.out.partitioned.collect(), ListPhenotypes.out.splitText( by: 100, file: true ), variant_reference_ch)

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
