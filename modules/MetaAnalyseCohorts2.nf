#!/bin/bash nextflow 


process MetaAnalyseCohorts2 {

    tag "chunk $chunk"

    publishDir "${params.outdir}", mode: 'copy'

    input:
      path genotype, stageAs: "genotypes_???"
   
    output:
      path 'MetaAnalysisResultsEncoded/*.parquet'

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    echo !{genotype}

    '''
}