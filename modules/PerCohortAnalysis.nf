#!/bin/bash nextflow 


process PerCohortAnalysis {

    tag {"chunk $chunk"}

    input:
      val chunk
      path mapper
      val th
      val nr_chunks
      path snp_inclusion
      path gene_inclusion
      path gene_per_cohort
      path covariate_filtering
      path genotype, stageAs: "genotypes_???"
      path expression, stageAs: "expression_???"
      path partial_derivatives, stageAs: "pd_???"
      val cohort
      val encoded

    output:
      tuple val(chunk), path('MetaAnalysisResultsEncoded/meta'), emit: meta
      tuple val(chunk), path('MetaAnalysisResultsEncoded/cohort'), emit: cohort

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    echo !{snp_inclusion}

    cohort=!{cohort.join(" ")}
    encoded=!{encoded.join(" ")}

    genotype=!{genotype.join(" ")}
    expression=!{expression.join(" ")}
    partial_derivatives=!{partial_derivatives.join(" ")}

    snp_inclusion=!{snp_inclusion.join(" ")}
    gene_inclusion=!{gene_inclusion.join(" ")}

    python2 -u !{baseDir}/bin/hase/hase.py \
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
    -ph_id_log !{gene_per_cohort} \
    -ci !{covariate_filtering}

    # Run combine command to combine the parquet files for every gene into one parquet file.
    # This allows using larger row group size, improving storage characteristics
    python2 !{baseDir}/bin/combine27.py \
    --path MetaAnalysisResultsEncoded/cohort \
    --phenotypes !{genes.join(" ")} \
    --cohorts ${cohort}
    '''
}