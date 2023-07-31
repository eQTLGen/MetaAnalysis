#!/bin/bash nextflow


process Partition {

    input:
      set val(chunk), path(parquet)

    output:
      path "partitioned_${chunk}", emit: partitioned
      file "phenotypes.txt", emit: phenotypes
      set val(chunk), val(1), emit: signal

    shell:
    '''
    mkdir "partitioned_!{chunk}"

    python2 -u !{baseDir}/bin/partition.py \
    --path !{parquet} \
    --out "partitioned_!{chunk}" \
    --out-list "phenotypes.txt"
    '''
}


process CleanPartition {
    tag {CleanPartition}

    input:
        set val(chunk), val(files_list)

    script:
    """
    clean_work_files.sh "${files_list[0]}"
    """
}