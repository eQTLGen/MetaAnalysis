#!/bin/bash nextflow


process Combine {
    cache true
    publishDir "${params.outdir}/eqtls", mode: 'copy', overwrite: true
    container 'quay.io/cawarmerdam/eqtlgen_phase2:latest'

    input:
      path partitioned
      val phenotype
      path reference

    output:
      tuple val(phenotype), path("phenotype*"), emit: parquet
      val 1, emit: signal

    shell:
    '''
    # Combining all parquet files for a specific phenotype
    # Allow multiple phenotypes to be processed in one go
    # Merge with reference for chrom to quickly be able to select variants in a specific chromosome

    python3 -u !{baseDir}/bin/combine.py \
    --path !{partitioned} \
    --pheno !{phenotype} \
    --out "." \
    --ref !{reference}
    '''
}

process CleanCombine {
    tag {CleanCombine}

    input:
        tuple val(phenotype), val(files_list)

    script:
    """
    clean_work_files.sh "${files_list[0]}"
    """
}
