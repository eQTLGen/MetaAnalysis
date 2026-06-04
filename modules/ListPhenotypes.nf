#!/bin/bash nextflow


process ListPhenotypes {
    cache true
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    executor 'local'

    input:
      path 'phenotype_list?.txt'

    output:
      path "phenotypes_unique.txt"

    shell:
    '''
    sort phenotype_list*.txt | uniq > "phenotypes_unique.txt"
    '''
}
