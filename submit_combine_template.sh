#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="MetaAnalysis"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load jdk/16.0.1
module load openjdk/11.0.2
module load any/singularity
module load squashfs/4.4

# Define paths
nextflow_path=/gpfs/space/GI/eQTLGen/EstBB_testing/MetaAnalysis/tools

set -f

input_path='/gpfs/space/GI/eQTLGen/freeze1/eqtl_mapping/output/empirical_4GenPC20ExpPC_2022-11-14/MetaAnalysisResultsEncoded/node_1_*_result.parquet'
output_folder=../output/MetaAnalysisResultsPartitioned

NXF_VER=21.10.6 ${nextflow_path}/nextflow run OutputPerPhenotype.nf \
--input ${input_path} \
--outdir ${output_folder} \
-resume \
-profile slurm,singularity
