#!/bin/bash nextflow


process PreMetaGeneFilter {
    publishDir "${params.outdir}/eqtls/", mode: 'copy', overwrite: true

    input:
        path gene_inclusion, stageAs: "gene_inclusion_???", arity: '1..*'
        val expression_platform
        path gene_custom
        path expression, stageAs: "expression_???", arity: '1..*'

    output:
        path "filtered_gene_list_2.txt"

    shell:
        // Get genes available for all included cohorts
        '''
        # 1
        python2 !{baseDir}/bin/pre_meta_gene_filter.py --expression !{expression.join(" ")} \
        --gene-inclusion !{gene_inclusion.join(" ")} --expression-platform '!{expression_platform.join("' '")}' \
        --output "filtered_gene_list.txt"

        cp "filtered_gene_list.txt" "filtered_gene_list_2.txt"

        if [ -f !{gene_custom} ]; then
          grep -f !{gene_custom} "filtered_gene_list.txt" > "filtered_gene_list_2.txt"
        fi
        '''
}
