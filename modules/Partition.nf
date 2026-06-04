#!/bin/bash nextflow


process Partition {
    container 'quay.io/cawarmerdam/eqtlgen_phase2:latest'

    input:
      tuple val(chunk), path(parquet)

    output:
      path "partitioned_${chunk}", emit: partitioned
      path "phenotypes.txt", emit: phenotypes
      tuple val(chunk), val(1), emit: signal

    shell:
    '''
    mkdir "partitioned_!{chunk}"

    python2 -u !{baseDir}/bin/partition.py \
    --path !{parquet} \
    --out "partitioned_!{chunk}" \
    --out-list "phenotypes.txt"
    '''
}


process OldPartitionPerCohort {
    cache true
    container 'quay.io/cawarmerdam/eqtlgen_phase2:latest'

    input:
      path parquet
      val chunk

    output:
      path "partitioned_${chunk}", emit: partitioned
      path "phenotypes.txt", emit: phenotypes

    shell:
    '''
    set +f

    mkdir "partitioned_!{chunk}"

    python3 -u !{baseDir}/bin/partition.py \
    --path !{parquet}/node_!{chunk}_*_per_cohort.parquet \
    --out "partitioned_!{chunk}" \
    --out-list "phenotypes.txt"
    '''
}


process CleanPartition {
    tag {CleanPartition}

    input:
        tuple val(chunk), val(files_list)

    script:
    """
    clean_work_files.sh "${files_list[0]}"
    """
}
