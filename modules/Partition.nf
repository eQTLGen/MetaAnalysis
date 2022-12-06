#!/bin/bash nextflow


process Partition {

    input:
      path parquet
      val chunk

    output:
      path "partitioned"

    shell:
    '''
    mkdir "partitioned"

    python2 -u !{baseDir}/bin/partition.py \
    --path !{parquet}/node_!{chunk}_*_result.parquet \
    --out "partitioned_!{chunk}"
    '''
}
