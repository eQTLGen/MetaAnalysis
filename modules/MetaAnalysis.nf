#!/bin/bash nextflow 


process PerCohortAnalysisPerGene {
    publishDir "${params.outdir}/eqtls/meta", mode: 'move', overwrite: true, pattern: 'MetaAnalysisResultsEncoded/meta/*', saveAs: { fn -> file(fn).getName() }
    publishDir "${params.outdir}/eqtls/cohort", mode: 'move', overwrite: true, pattern: 'MetaAnalysisResultsEncoded/cohort/*', saveAs: { fn -> file(fn).getName() }
    scratch true

    input:
      val th
      path genes
      path variants_per_cohort
      path mapper
      path covariate_filtering
      val cohort
      val encoded
      path genotype, stageAs: "genotypes_???", arity: '1..*'
      path expression, stageAs: "expression_???", arity: '1..*'
      path partial_derivatives, stageAs: "pd_???", arity: '1..*'
      path snp_inclusion, stageAs: "snp_inclusion_???", arity: '1..*'
      path gene_inclusion, stageAs: "gene_inclusion_???", arity: '1..*'

    output:
      tuple path('MetaAnalysisResultsEncoded/meta/*'), emit: meta
      tuple path('MetaAnalysisResultsEncoded/cohort/*'), emit: cohort

    shell:
    snps_per_cohort_arg = (variants_per_cohort.name != 'NO_FILE') ? "-snp_id_log ${variants_per_cohort}" : ""

    snp_inclusion_per_cohort = snp_inclusion.join(' ')
    if (variants_per_cohort.name != 'NO_FILE') {
        snp_inclusion_per_cohort = snp_inclusion.name.collect { filename -> "intersect_$filename" }.join(' ')
    }
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    mkdir tmp_files

    cp -r genotypes* tmp_files/
    cp -r expression* tmp_files/
    cp -r pd* tmp_files/

    if [[ !{variants_per_cohort.name} != 'NO_FILE' ]]; then
        # Filter snp inclusion files to only contain snps to be included
        for snp_inclusion_file in !{snp_inclusion.join(' ')}; do
          echo "ID" > "intersect_${snp_inclusion_file}"
          comm -12 <(tail -n +2 !{variants_per_cohort} | sort) <(sort ${snp_inclusion_file}) >> "intersect_${snp_inclusion_file}"
        done
    fi

    # Filter gene inclusion files to only contain genes to be ran in this chunk
    for gene_inclusion_file in !{gene_inclusion.join(' ')}; do
      echo "ID" > "intersect_${gene_inclusion_file}"
      comm -12 <(tail -n +2 !{genes} | sort) <(sort ${gene_inclusion_file}) >> "intersect_${gene_inclusion_file}"
    done

    python2 -u !{baseDir}/bin/hase/hase.py \
      -study_name !{cohort.join(" ")} \
      -g !{genotype.name.collect{filename -> "tmp_files/$filename"}.join(' ')} \
      -ph !{expression.name.collect{filename -> "tmp_files/$filename"}.join(' ')}  \
      -derivatives !{partial_derivatives.name.collect{filename -> "tmp_files/$filename"}.join(' ')} \
      -mapper !{mapper}/ \
      -o MetaAnalysisResultsEncodedTmp \
      -mode meta-classic \
      -encoded !{encoded.join(" ")} \
      -max-missingness-rate 0.8 \
      -thr !{th} \
      -mapper_chunk 500 \
      -ref_name 1000G-30x_ref \
      -snp_id_inc !{snp_inclusion_per_cohort} \
      -ph_id_inc !{gene_inclusion.name.collect { filename -> "intersect_$filename" }.join(' ')} \
      !{snps_per_cohort_arg} \
      -ci !{covariate_filtering}

    # Run combine command to combine the parquet files for every gene into one parquet file.
    # This allows using larger row group size, improving storage characteristics
    python2 !{baseDir}/bin/combine27.py \
    --path MetaAnalysisResultsEncodedTmp/cohort \
    --out MetaAnalysisResultsEncoded/cohort \
    --phenotypes !{genes} \
    --cohorts !{cohort.join(" ")}

    python2 !{baseDir}/bin/combine27.py \
    --path MetaAnalysisResultsEncodedTmp/meta \
    --out MetaAnalysisResultsEncoded/meta \
    --phenotypes !{genes}

    '''
}
