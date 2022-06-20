#!/usr/bin/env nextflow

def helpmessage() {

log.info"""

HASE meta-analyzer ~ v${workflow.manifest.version}"
=================================================
Pipeline for running encoded meta-analysis over multiple studies in the central site.

This pipeline mostly recapitulates the steps as outlined in HASE wiki page (https://github.com/roshchupkin/hase/wiki/HD-Meta-analysis-GWAS), with some bugfixes. Pay attention that it does not explicitly overwrite the output folder, so clean it before re-running the pipeline.

Usage:

nextflow run HaseMetaAnalysis.nf \
--inputpath [path to organised inputs] \
--mapperpath [path to mapper files]\
--outputpath [path to results folder] \

Mandatory arguments:
--inputpath       Path to organised encoded data.
--mapperpath      Path to the mappers folder, containing the mapper files for each cohort.
--outputpath      Path to results folder where meta-analysis outputs are written.

Optional arguments:
--maf            MAF threshold for SNPs to include into combined analysis. Defaults to 0.
--th             T-statistic threshold for writing out the results. Defaults to 0.
--snps           SNPs to include to the analysis. Files with header ID and SNP IDs on each row, or file with header CHR,POS and SNP positions in format [1,1165539] on each row. To use it, [HASE path/data] needs to include *.ref_info.h5 file.
--phenotypes     Phenotypes to include to the analysis. Files with header ID and phenotype names on each row.
--chunks         In how many chunks to run the analysis. Defaults to 100.
--covariates     What covariates to include into the analysis for each cohort. If specified, the it has to be the file which has each line formatted as following: [cohort name] [index of covariate to include 1] [index of covariate to include 2] etc. Indexing is 1-based and defaults to all covariates.
--single         Whether the input contains data for only one cohort. yes/no, defaults to no.
--permuted       Whether the analysis should be done on permuted data. yes/no, defaults to no.

""".stripIndent()

}

//Default parameters
params.inputpath = ''
params.outputpath = ''
params.mapperpath = ''
params.maf = 0
params.th = 0
params.snps = ''
params.phenotypes = ''
params.covariates = ''
params.chunks = 100
params.single = 'no'
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
summary['Input folder']                             = params.inputpath
summary['Input mapper path']                        = params.mapperpath
summary['Output directory']                         = params.outputpath
summary['MAF threshold']                            = params.maf
summary['T-statistic threshold']                    = params.th
summary['SNP filter']                               = params.snps
summary['Phenotype filter']                         = params.phenotypes
summary['Covariate filter']                         = params.covariates
summary['Nr. of chunks']                            = params.chunks
summary['Single cohort']                            = params.single
summary['Permuted data']                            = params.permuted

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "================================================="

// Count the number of input datasets
//StudyNames = Channel.value(params.studynames).val.split('[ ]{1,}')
//NrOfDatasets = StudyNames.size().value

// Input file paths
inputpath = file(params.inputpath)
mapper = file(params.mapperpath)

// Optional arguments
maf = Channel.value(params.maf)
th = params.th
chunks = params.chunks

if (params.phenotypes) { 
  phen_file = file(params.phenotypes) 
  phen_filtering = Channel.value("-ph_id_inc ${phen_file}")
} else {phen_filtering = Channel.value("")}

if (params.snps) { 
  snp_file = file(params.snps) 
  snp_filtering = Channel.value("-snp_id_inc ${snp_file}")
} else {snp_filtering = Channel.value("")}

if (params.covariates) { 
  covariate_file = file(params.covariates) 
  covariate_filtering = Channel.value("-ci ${covariate_file}")
} else {covariate_filtering = Channel.value("")}

process GetInputVariables {

  input:
    path(inputpath) from inputpath

  output:
    env files_paths
    env genotype_paths
    env phenotype_paths
    env pd_paths
    env cohort_names
    env encodedstatus

  shell:
  if (params.permuted == 'no')
      '''
    files_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/permuted/)
    genotype_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/permuted/EncodedGenotypeData/)
    phenotype_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/permuted/EncodedPhenotypeData/)
    pd_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/permuted/pd_shared/)

    cohort_names=$(sed 's/_IntermediateFiles.*//g' <<<"$files_paths") 
    cohort_names=$(sed 's/.*\\///g' <<<"$cohort_names")

    encodedstatus=$(sed 's/.*permuted\\//permuted/g' <<<"$files_paths")
    encodedstatus=$(sed 's/.*permuted/1/g' <<<"$encodedstatus")
    '''
  else
    '''
    files_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/empirical/)
    genotype_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/empirical/EncodedGenotypeData/)
    phenotype_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/empirical/EncodedPhenotypeData/)
    pd_paths=$(ls -d !{inputpath}/*/*_IntermediateFilesEncoded_to_upload/empirical/pd_shared/)

    cohort_names=$(sed 's/_IntermediateFiles.*//g' <<<"$files_paths") 
    cohort_names=$(sed 's/.*\\///g' <<<"$cohort_names")

    encodedstatus=$(sed 's/.*empirical\\//empirical/g' <<<"$files_paths")
    encodedstatus=$(sed 's/.*empirical/1/g' <<<"$encodedstatus")
    '''

}

chunk = Channel.from(1..chunks)

process MetaAnalyseCohorts {

    tag "chunk $chunk"

    publishDir "${params.outputpath}", mode: 'copy'

    input:
      val chunk from chunk
      path inputpath from inputpath
      path mapper from mapper
      val(maf) from maf
      val(th) from th
      val(nr_chunks) from chunks
      val(snp_filtering) from snp_filtering
      val(phen_filtering) from phen_filtering
      val(covariate_filtering) from covariate_filtering
      val(single) from params.single
      env files_paths
      env genotype_paths
      env phenotype_paths
      env pd_paths
      env cohort_names
      env encodedstatus

    output:
      path 'MetaAnalysisResultsEncoded/*.parquet' into MetaAnalysisResultsEncoded

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    python2 !{baseDir}/bin/hase/hase.py \
    -study_name ${cohort_names} \
    -g ${genotype_paths} \
    -ph ${phenotype_paths} \
    -derivatives ${pd_paths} \
    -mapper !{mapper}/ \
    -o MetaAnalysisResultsEncoded \
    -mode meta-classic \
    -encoded ${encodedstatus} \
    -allow_missingness \
    -node !{nr_chunks} !{chunk} \
    -maf !{maf} \
    -th !{th} \
    -mapper_chunk 500 \
    -ref_name 1000G-30x_ref \
    -cluster "y" \
    !{snp_filtering} \
    !{phen_filtering} \
    !{covariate_filtering}
    '''
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
