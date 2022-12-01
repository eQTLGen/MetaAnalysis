#!/bin/bash nextflow 


process MetaAnalyseCohorts {

    tag "chunk $chunk"

    publishDir "${params.outdir}", mode: 'copy'

    input:
      val chunk
      path mapper
      val th
      val nr_chunks
      path snp_inclusion
      path gene_inclusion
      path covariate_filtering 
      path genotype, stageAs: "genotypes_???"
      path expression, stageAs: "expression_???"
      path partial_derivatives, stageAs: "pd_???"
      val cohort
      val encoded

    output:
      path 'MetaAnalysisResultsEncoded/*.parquet'

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    echo !{snp_inclusion} 

    cohort=$(echo !{cohort} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')
    encoded=$(echo !{encoded} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')

    genotype=$(echo !{genotype} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')
    expression=$(echo !{expression} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')
    partial_derivatives=$(echo !{partial_derivatives} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')

    snp_inclusion=$(echo !{snp_inclusion} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')
    gene_inclusion=$(echo !{gene_inclusion} | sed -r 's/\\]//g' | sed -r 's/\\[//g' | sed -r 's/,//g')

    python2 !{baseDir}/bin/hase/hase.py \
    -study_name ${cohort} \
    -g ${genotype} \
    -ph ${expression} \
    -derivatives ${partial_derivatives} \
    -mapper !{mapper}/ \
    -o MetaAnalysisResultsEncoded \
    -mode meta-classic \
    -encoded ${encoded} \
    -allow_missingness \
    -node !{nr_chunks} !{chunk} \
    -th !{th} \
    -mapper_chunk 500 \
    -ref_name 1000G-30x_ref \
    -cluster "y" \
    -snp_id_inc ${snp_inclusion} \
    -ph_id_inc ${gene_inclusion} \
    -ci !{covariate_filtering}
    '''
}
