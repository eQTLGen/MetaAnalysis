#!/bin/bash nextflow 

process SubsetGenesInclusion {

    input:
      path gene_inclusion
      path gene_per_cohort

    output:
      path 'intersect_*'

    shell:
    '''
    for gene_inclusion_file in !{gene_inclusion}; do
      echo "ID" > "intersect_${gene_inclusion_file}"
      comm -12 <(tail -n +2 !{gene_per_cohort} | sort) <(sort ${gene_inclusion_file}) >> "intersect_${gene_inclusion_file}"
    done
    '''
}
