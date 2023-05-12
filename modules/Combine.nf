#!/bin/bash nextflow


process Combine {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path partitioned
      val phenotype

    output:
      path "phenotype*"

    shell:
    '''as
    # Combining all parquet files for a specific phenotype

    python2 -u !{baseDir}/bin/combine.py \
    --path !{partitioned} \
    --pheno !{phenotype} \
    --out "." \
    --ref "!{baseDir}/bin/hase/data/1000G-30x.ref_info.h5"
    '''
}
