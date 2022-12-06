#!/bin/bash nextflow


process Partition {

    input:
      path parquet

    output:
      path "partitioned"

    shell:
    '''
    mkdir "partitioned"

    python2 !{baseDir}/bin/partition.py \
    --path !{parquet} \
    --out "partitioned"
    '''
}
