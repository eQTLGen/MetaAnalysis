#!/bin/bash nextflow


process PreMetaGeneFilter {
    publishDir "${params.outdir}/eqtls/", mode: 'copy', overwrite: true

    input:
        path gene_inclusion
        path gene_custom
        path expression

    output:
        path "filtered_gene_list.txt"

    shell:
        // Get genes available for all included cohorts
        '''
        pre_meta_gene_filter.py --expression !{expression.join(" ")} --gene-inclusion !{(gene_custom.name != 'NO_FILE') ? gene_custom} !{gene_inclusion.join(" ")} \
        --output "filtered_gene_list.txt"
        '''