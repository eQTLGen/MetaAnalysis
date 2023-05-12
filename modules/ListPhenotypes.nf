#!/bin/bash nextflow


process ListPhenotypes {

    input:
      path partitioned

    output:
      path "phenotypes.txt"

    shell:
    '''
    ls !{partitioned}/phenotype_* | basename | awk 'fs="_" { print $2 }' > "phenotypes.txt"
    '''
}
