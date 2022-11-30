#!/bin/bash nextflow 


process Combine {

    input:
      path partitioned
      val phenotype
      path out

    shell:
    '''
    python2 !{baseDir}/bin/combine.py \
    --path ${parquet} \
    --pheno ${phenotype} \
    --out ${out} \
    '''
}
