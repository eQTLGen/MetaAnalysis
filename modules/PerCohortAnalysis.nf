#!/bin/bash nextflow 


process PerCohortAnalysis {

    tag {"chunk $chunk"}

    input:
      val th
      val genes
      path mapper
      path covariate_filtering
      tuple val(cohort),
            val(encoded),
            path(genotype, stageAs: "genotypes_???"),
            path(expression, stageAs: "expression_???"),
            path(partial_derivatives, stageAs: "pd_???"),
            path(snp_inclusion, stageAs: "snp_inclusion_???"),
            path(gene_inclusion, stageAs: "gene_inclusion_???")

    output:
      tuple val(chunk), path('MetaAnalysisResultsEncoded/meta'), emit: meta
      tuple val(chunk), path('MetaAnalysisResultsEncoded/cohort'), emit: cohort

    shell:
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    echo !{snp_inclusion}

    echo !{genes.join("\n")} > gene_chunk.txt

    for gene_inclusion_file in !{gene_inclusion.join(" ")}; do
      echo "ID" > "intersect_${gene_inclusion_file}"
      comm -12 gene_chunk.txt | sort) <(sort ${gene_inclusion_file}) >> "intersect_${gene_inclusion_file}"
    done

    cohort=!{cohort.join(" ")}
    encoded=!{encoded.join(" ")}

    genotype=!{genotype.name.join(" ")}
    expression=!{expression.name.join(" ")}
    partial_derivatives=!{partial_derivatives.name.join(" ")}

    snp_inclusion=!{snp_inclusion.name.join(" ")}
    gene_inclusion=!{gene_inclusion.name.collect { filename -> "intersect_$filename" }.join(" ")}

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
    -ph_id_log ${gene_inclusion} \
    -ci !{covariate_filtering}

    # Run combine command to combine the parquet files for every gene into one parquet file.
    # This allows using larger row group size, improving storage characteristics
    python2 !{baseDir}/bin/combine27.py \
    --path MetaAnalysisResultsEncoded/cohort \
    --phenotypes !{genes.join(" ")} \
    --cohorts ${cohort}
    '''
}