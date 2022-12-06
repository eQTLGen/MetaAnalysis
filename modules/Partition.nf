#!/bin/bash nextflow


process Partition {

    input:
      path parquet
      val chunk

    output:
      path "partitioned_!{chunk}"

    shell:
    '''
    mkdir "partitioned_!{chunk}"

    python2 -u !{baseDir}/bin/partition.py \
    --path !{parquet}/node_!{chunk}_*_result.parquet \
    --out "partitioned_!{chunk}"
    '''
}
