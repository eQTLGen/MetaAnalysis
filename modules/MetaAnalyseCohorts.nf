#!/bin/bash nextflow


process MetaAnalyseCohortsPerGene {

    input:
      val th
      val genes
      path mapper
      path covariate_filtering
      tuple val(cohort),
            val(encoded),
            path(genotype, stageAs: "genotypes_???"),
            path(expression, stageAs: "expression_???"),
            path(partial_derivatives, stageAs: "pd_???")
            path(snp_inclusion, stageAs: "snp_inclusion_???"),
            path(gene_inclusion, stageAs: "gene_inclusion_???"),

    output:
      tuple path('MetaAnalysisResultsEncoded/meta'), emit: meta

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

    genotype=!{genotype.join(" ")}
    expression=!{expression.join(" ")}
    partial_derivatives=!{partial_derivatives.join(" ")}

    snp_inclusion=!{snp_inclusion.join(" ")}
    gene_inclusion=!{gene_inclusion.collect { "intersect_$gene_inclusion_file" }.join(" ")}

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
    -th !{th} \
    -mapper_chunk 500 \
    -ref_name 1000G-30x_ref \
    -snp_id_inc ${snp_inclusion} \
    -ph_id_inc ${gene_inclusion} \
    -ci !{covariate_filtering}
    '''
}

process CleanMetaAnalyseCohorts {
    tag {CleanMetaAnalyseCohorts}

    input:
        tuple val(id), val(files_list)

    script:
    """
    #clean_work_files.sh "${files_list[0]}"
    """
}

