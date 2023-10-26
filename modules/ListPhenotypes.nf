#!/bin/bash nextflow


process ListPhenotypes {

    input:
      path partitioned

    output:
      path "phenotypes.txt"

    shell:
    '''
    sort phenotype_list*.txt | uniq > "phenotypes_unique.txt"
    '''
}
