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
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# Define paths 
nextflow_path=[folder where Nextflow executable is]

input_path=[path to input files]
output_folder=[path to the folder where output files are written]

NXF_VER=20.10.6 ${nextflow_path}/nextflow run OutputPerPhenotype.nf \
--input ${input_path} \
--outdir ${output_folder} \
-resume \
-profile slurm,singularity
