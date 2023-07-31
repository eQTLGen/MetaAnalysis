#!/bin/bash nextflow


process Combine {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path partitioned
      val phenotype

    output:
      set val(phenotype), path("phenotype*"), emit: parquet
      val 1, emit: signal

    shell:
    '''
    # Combining all parquet files for a specific phenotype

    python2 -u !{baseDir}/bin/combine.py \
    --path !{partitioned} \
    --pheno !{phenotype} \
    --out "." \
    --ref "!{baseDir}/bin/hase/data/1000G-30x.ref_info.h5"
    '''
}

process CleanCombine {
    tag {CleanCombine}

    input:
        set val(phenotype), val(files_list)

    script:
    """
    clean_work_files.sh "${files_list[0]}"
    """
}