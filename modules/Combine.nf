#!/bin/bash nextflow 


process Combine {

    input:
      path partitioned
      val phenotype
      path out

    shell:
    '''
    python2 -u !{baseDir}/bin/combine.py \
    --path !{partitioned} \
    --pheno !{phenotype} \
    --out !{out} \
    '''
}
