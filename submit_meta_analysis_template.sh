#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="RunHASEwNextflow"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load java-1.8.0_40
module load python/2.7.15/native
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=[Path to tool folder with nextflow]

NXF_VER=20.10.0 ${nextflow_path}/nextflow run HaseMetaAnalysis.nf \
--genopath '[Path(s) to encoded genotype folder(s)]' \
--expressionpath '[Path(s) to encoded expression data folder(s)]' \
--mapperpath '[Path to folder with mapper files from all datasets]' \
--pdpath '[Path(s) to partial derivatives folder(s)]' \
--studynames '[CohortName_GeneExpressionPlatform_for_first_dataset CohortName_GeneExpressionPlatform_for_second_dataset]' \
--outputpath '[Path to results output folder]' \
--maf 0.01 \
-with-report Meta_analysis_report.html \
-resume \
-profile singularity_profile,cluster_slurm
