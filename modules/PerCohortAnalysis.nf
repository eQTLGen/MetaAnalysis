#!/bin/bash nextflow 


process PerCohortAnalysis {

    tag {"chunk $chunk"}

    input:
      val th
      path genes
      path mapper
      path covariate_filtering
      val cohort
      val encoded
      path genotype, stageAs: "genotypes_???"
      path expression, stageAs: "expression_???"
      path partial_derivatives, stageAs: "pd_???"
      path snp_inclusion, stageAs: "snp_inclusion_???"
      path gene_inclusion, stageAs: "gene_inclusion_???"

    output:
      tuple val(chunk), path('MetaAnalysisResultsEncoded/meta'), emit: meta
      tuple val(chunk), path('MetaAnalysisResultsEncoded/cohort'), emit: cohort

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    for gene_inclusion_file in !{gene_inclusion.join(' ')}; do
      echo "ID" > "intersect_${gene_inclusion_file}"
      sort ${gene_inclusion_file} > sorted_gene_inclusion.txt
      comm -12 <(tail -n +2 !{genes} | sort) <(sort ${gene_inclusion_file}) >> "intersect_${gene_inclusion_file}"
    done

    python2 -u !{baseDir}/bin/hase/hase.py \
      -study_name !{cohort.join(" ")} \
      -g !{genotype.name.join(" ")} \
      -ph !{expression.name.join(" ")} \
      -derivatives !{partial_derivatives.name.join(" ")} \
      -mapper !{mapper}/ \
      -o MetaAnalysisResultsEncoded \
      -mode meta-classic \
      -encoded !{encoded.join(" ")} \
      -allow_missingness \
      -th !{th} \
      -mapper_chunk 500 \
      -ref_name 1000G-30x_ref \
      -snp_id_inc !{snp_inclusion.name.join(" ")} \
      -ph_id_inc !{gene_inclusion.name.collect { filename -> "intersect_$filename" }.join(' ')} \
      -ph_id_log !{gene_inclusion.name.collect { filename -> "intersect_$filename" }.join(' ')} \
      -ci !{covariate_filtering}

    # Run combine command to combine the parquet files for every gene into one parquet file.
    # This allows using larger row group size, improving storage characteristics
    python2 !{baseDir}/bin/combine27.py \
    --path MetaAnalysisResultsEncoded/cohort \
    --phenotypes !{genes} \
    --cohorts ${cohort}
    '''
}
