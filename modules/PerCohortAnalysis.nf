#!/bin/bash nextflow 


process PerCohortAnalysis {

    tag {"chunk $chunk"}

    input:
      val th
      path genes
      path genes_per_cohort
      path variants_per_cohort
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
    snps_per_cohort_arg = (variants_per_cohort.name != 'NO_FILE') ? "-snp_id_log ${variants_per_cohort}" : ""

    snp_inclusion_per_cohort = snp_inclusion.join(' ')
    if (variants_per_cohort.name != 'NO_FILE') {
        snp_inclusion_per_cohort = snp_inclusion.name.collect { filename -> "intersect_$filename" }.join(' ')
    }
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    # Filter gene inclusion files to only contain genes to be ran in this chunk
    for gene_inclusion_file in !{gene_inclusion.join(' ')}; do
      echo "ID" > "intersect_${gene_inclusion_file}"
      comm -12 <(tail -n +2 !{genes} | sort) <(sort ${gene_inclusion_file}) >> "intersect_${gene_inclusion_file}"
    done

    # Filter snp inclusion files to only contain snps to be included
    for snp_inclusion_file in !{snp_inclusion.join(' ')}; do
      echo "ID" > "intersect_${snp_inclusion_file}"
      comm -12 <(tail -n +2 !{variants_per_cohort} | sort) <(sort ${snp_inclusion_file}) >> "intersect_${snp_inclusion_file}"
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
      -snp_id_inc !{snp_inclusion_per_cohort} \
      -ph_id_inc !{gene_inclusion.name.collect { filename -> "intersect_$filename" }.join(' ')} \
      -ph_id_log !{genes_per_cohort} \
      !{snps_per_cohort_arg} \
      -ci !{covariate_filtering}

    # Run combine command to combine the parquet files for every gene into one parquet file.
    # This allows using larger row group size, improving storage characteristics
    python2 !{baseDir}/bin/combine27.py \
    --path MetaAnalysisResultsEncoded/cohort \
    --phenotypes !{genes} \
    --cohorts ${cohort}
    '''
}
