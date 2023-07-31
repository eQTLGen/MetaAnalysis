#!/bin/bash nextflow


process PreparePermuted {

    publishDir "${params.outputpath}", pattern: "merged_genotypes_${cohort}", mode: "copyNoFollow"

    input:
      tuple val(cohort), path(empirical, stageAs: 'empirical'), path(permuted, stageAs: 'permuted')

    output:
      path "merged_genotypes_${cohort}"

    shell:
    '''
      out="merged_genotypes_!{cohort}"
      mkdir $out

      ln -s $(readlink -f "!{permuted}/individuals") "${out}/individuals"
      ln -s $(readlink -f "!{permuted}/genotype") "${out}/genotype"

      ln -s $(readlink -f "!{empirical}/probes") "${out}/probes"
    '''
}
