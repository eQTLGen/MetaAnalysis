#!/bin/bash nextflow


process Combine {

    publishDir "${params.outdir}", mode: 'move', overwrite: false

    input:
      path partitioned
      val phenotype

    output:
      path "output/*"

    shell:
    '''
    python2 -u !{baseDir}/bin/combine.py \
    --path !{partitioned} \
    --pheno !{phenotype} \
    --out "output" \
    '''
}
