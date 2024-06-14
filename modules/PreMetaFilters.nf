#!/bin/bash nextflow


process PreMetaGeneFilter {
    publishDir "${params.outdir}/eqtls/", mode: 'copy', overwrite: true

    input:
        path gene_inclusion, stageAs: "gene_inclusion_???", arity: '1..*'
        path gene_custom
        path expression, stageAs: "expression_???", arity: '1..*'

    output:
        path "filtered_gene_list_2.txt"

    shell:
        // Get genes available for all included cohorts
        '''
        python2 !{baseDir}/bin/pre_meta_gene_filter.py --expression !{expression.join(" ")} --gene-inclusion !{gene_inclusion.join(" ")} \
        --output "filtered_gene_list.txt"

        grep -f !{gene_custom} "filtered_gene_list.txt" > "filtered_gene_list_2.txt"
        '''
}
