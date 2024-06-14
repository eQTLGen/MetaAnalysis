#!/bin/bash nextflow


process Combine {
    containerOptions "--bind ${outdir}:${outdir}"

    input:
      path genes
      val outdir

    shell:
    '''
    # Combining all parquet files for a specific phenotype
    # Allow multiple phenotypes to be processed in one go

    python2 !{baseDir}/bin/combine27.py \
        --path !{outdir}/meta \
        --out-dir !{outdir}/meta \
        --out-tag combined \
        --from-parquet \
        --phenotypes !{genes}
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
