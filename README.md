## Encoded meta-analysis pipeline

This is Nextflow pipeline for running encoded QTL meta-analysis on the encoded data matrices.

### Usage instructions

#### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=8 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

#### Setup of the pipeline
You can either clone it by using git (if available in HPC):

`git clone https://github.com/eQTLGen/MetaAnalysis.git`

Or just download this from the gitlab/github download link and unzip.

#### Input files

- It is expected that all encoded input datasets are prepared with eQTLGen phase II pipelines, as explained in [eQTLGen phase 2 cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html).

- Paths to the input files are collected to tab-separated **master table**. Each row in the master table corresponds to one dataset in the encoded meta-analysis. You can initially prepare this in Excel and it must have the following structure and column names:
  - cohort: Name of the cohort.
  - genotype: Path to encoded genotype folder.
  - expression: Path to encoded phenotype folder.
  - partial_derivatives: Path to partial derivatives folder.
  - encoded: Indicator whether dataset is encoded, must be 1 for encoded dataset.
  - snp_inclusion: File with SNPs which pass QC filters and should be included to the analysis for this cohort. See below how to make it.
  - gene_inclusion: File with genes which pass QC filters and should be included to the analysis for this cohort. See below how to make it.

- Separate **covariate file**: this file specifies which covariates from partial derivative data should be included to the analysis for each cohort. See below how to make it.

- Folder with **mapper files**: this must include mapper files for all cohorts included to the analysis. Just copy the mapper files from the prepared data of each included cohort into one folder.
  
#### Usage instructions

##### Initial preparations

We provide separate R script which can be used to easily make per-cohort SNP inlcusion lists, gene inclusion lists and covariate file. The script template is findable from `bin/make_inclusion_lists.R`. This script expects that data is prepared as specified in [eQTLGen phase 2 cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html).

You need following R packages to run this helper script: `data.table`, `stringr`, `fs`, `tidyr`, `gtools`.

Standard QC applies following filters, you can also adjust the settings:
- SNP QC:
  - MAF>0.01
  - Hardy-Weinberg disequilibrium P>1e-6
  - Imputation Mach R2>0.4
  - Call rate=1 (for imputed data this is always the case)
- Gene expression QC:
  - Per gene require that at least 80% of normalised and log2-transformed values are unique. This ensures that gene is expressed (majority of values are not 0) and shows at least some variability.
- Covariates:
  - By default add 4 first genotype PCs as covariates (correct for population stratification).
  - By default add 20 first expression PCs as covariates (correct for non-genetic expression variance).
  - When setting `gen_reg_pcs <- TRUE`, it does not add those expression PCs which are under genetic control. If this is active, script checks whether any of the first 20 expression PCs are under genetic control and includes only those which are not.
  
  Optimal covariate settings depend on the properties of the datasets.

```r
library(data.table)
library(stringr)
library(fs)
library(tidyr)
library(gtools)

setDTthreads(8)

#########################
# arguments and settings#
#########################

outputpath <- "[path to output folder where SNP and gene inclusion files are written]"

dataset_path <- c("/[path to cohort 1]/[name of cohort 1]_IntermediateFilesEncoded_to_upload/",
"/[path to cohort 2]/[name of cohort 2]_IntermediateFilesEncoded_to_upload/",
"/[path to cohort 3]/[name of cohort 3]_IntermediateFilesEncoded_to_upload/")

# You can adjust the SNP QC thresholds here
MAF_thresh <- 0.01
HWE_thresh <- 1e-6
MACH_R2_thresh <- 0.4
CALL_thresh <- 1

# You can adjust the gene QC thresholds here
gene_uniq_thresh <- 0.8

# You can adjust the P-value threshold for genetically regulated exp PCs here
GenRegPC_thresh <- 1e-7

# Number of first genetic PCs to add as covariates
gen_pcs <- 4
# Number of first expression PCs to add as covariates
exp_pcs <- 20
# Shall we remove those expression PCs which are under genetic regulation? By default not.
gen_reg_pcs <- FALSE

####################################

# Make SNP QC files
fileConn<-file(paste0(outputpath, "/covariate_indices.txt"), "w")

sink(paste0(outputpath, "/filter_logs.log"))
cat("Dataset\tN\tNr. SNPs\tFiltered SNPs\tTyped SNPs\tSNPs AC>5; not filtered\tGenes original\tGenes filtered\tGenetically regulated PCs\n")

for (dataset in dataset_path){

dataset_name <- str_replace(as.character(fs::dir_ls(path = paste0(dataset, "/SumStats/genotypes/SNPQC/"), glob = "*_SNPQC.txt.gz")), ".*/SNPQC/", "")
dataset_name <- str_replace(dataset_name, "_SNPQC.txt.gz", "")

message(paste0("Processing: ", dataset_name))
message("Make SNP QC files")

snpqc_filename <- fs::dir_ls(path = paste0(dataset, "/SumStats/genotypes/SNPQC/"), glob = "*_SNPQC.txt.gz")
snpqc <- fread(snpqc_filename, showProgress = FALSE)

# filters
orig_snps <- nrow(snpqc)
snpqc <- snpqc[MAF > MAF_thresh & HWE > HWE_thresh & MACH_R2 > MACH_R2_thresh & CALL == CALL_thresh]

# Show variants which have at least 5 allele counts for rarer genotype
hom1_gen_counts <- str_replace(snpqc$Genotype_Counts, ", .*", "")
snpqc$hom1_gen_counts <- as.numeric(str_replace(hom1_gen_counts, ".*: ", ""))

hom2_gen_counts <- str_replace(snpqc$Genotype_Counts, ".*: ", "")
snpqc$hom2_gen_counts <- as.numeric(str_replace(hom2_gen_counts, ",", ""))

# Calculate N
N_split <- unlist(str_split(snpqc$Genotype_Counts[1], ", "))
N_split <- str_replace(N_split, ".*: ", "")
N_split <- as.numeric(str_replace(N_split, ",", ""))
N_split <- sum(N_split)

fwrite(data.table(ID = snpqc$ID), paste0(outputpath, "/", dataset_name, "_SnpsToInclude.txt"), sep = "\t", row.names = FALSE)

geneqc_processed <- fread(paste0(dataset, "/SumStats/expression/processed_gene_summary.txt"))
genes_to_include <- geneqc_processed[geneqc_processed$nr_unique_values > gene_uniq_thresh * max(geneqc_processed$nr_values)]$gene
fwrite(data.frame(ID = genes_to_include), paste0(outputpath, "/", dataset_name, "_GenesToInclude.txt"), sep = "\t", row.names = FALSE)

message("Check genetically regulated expression PCs")
# Report PCs under genetic regulation, make covariate indices file
genetic_pcs <- fread(paste0(dataset, "/", dataset_name, "_GenRegPcs.txt.gz"))
genetic_pcs <- genetic_pcs[SNP %in% snpqc$ID & `P-value` < GenRegPC_thresh]
genetic_pcs <- mixedsort(unique(genetic_pcs$Phenotype))

# Identify which covariates to add
help_table <- data.table(covariate = c(paste0("GenPC", 1:10), paste0("ExpPC", 1:100)), index = 1:110)
# Filter in the covariates to include
include_table_gen <- head(help_table[str_detect(help_table$covariate, "GenPC"), ], gen_pcs)
include_table_exp <- head(help_table[str_detect(help_table$covariate, "ExpPC"), ], exp_pcs)
include_table <- rbind(include_table_gen, include_table_exp)

if(gen_reg_pcs == TRUE){
    include_table <- include_table[!include_table$covariate %in% genetic_pcs, ]
    }

write(paste0(c(dataset_name, include_table$index), collapse = "\t"), file = fileConn, append = TRUE)

# Write log
cat(paste(dataset_name, 
N_split, 
orig_snps,
nrow(snpqc), 
nrow(snpqc[snpqc$TYPED == "1",]), 
nrow(snpqc[hom1_gen_counts > 5 & hom2_gen_counts > 5]),
nrow(geneqc_processed),
length(genes_to_include),
paste(genetic_pcs, collapse = ", "),
"\n", sep = "\t"))

}

message("Finalising...")
sink()

```

The output of this script writes into specified folder for each dataset:
- SNP inclusion files (`*_SnpsToInclude.txt`)
- Gene inclusion files (`*_GenesToInclude.txt`)

It also writes covariate file to the same folder: `covariate_indices.txt`

Additionally it writes log file which gives some summary statistics about each cohort: `filter_logs.log`. This one can be used to identify datasets which behave differently from the rest.

#### Running the meta-analysis command

Go to folder `MetaAnalysis` and modify the Slurm script template `submit_meta_analysis_template.sh` with your input paths. Below is an example template for Slurm scheduler.

```bash
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
nextflow_path=[folder where Nextflow executable is] # folder where Nextflow executable is

mastertable=[path to mastertable]
mapper_folder=[path to the folder with all mapper files]
output_folder=[path to the folder where output files are written]


NXF_VER=20.10.6 ${nextflow_path}/nextflow run HaseMetaAnalysis.nf \
--mastertable ${mastertable} \
--mapperpath ${mapper_folder} \
--covariates /gpfs/space/GI/eQTLGen/EstBB_testing/MetaAnalysis/helpfiles2/covariate_indices.txt \
--chunks 100 \
--outdir ${output_folder} \
-resume \
-profile slurm,singularity
```

You can save the modified script version with different name, e.g. `submit_meta_analysis.sh`.

Then submit the job `sbatch submit_meta_analysis.sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate work directory is made to the folder and contains all interim files.

##### Monitoring and debugging
- Monitoring:
  - Monitor the slurm-***.out log file and check if all the steps finish without error. Trick: command watch tail -n 20 slurm-***.out helps you to interactively monitor the status of the jobs.
  - Use squeue -u [YourUserName] to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file is automatically written to your output folder pipeline_info subfolder, for potential errors or warnings. E.g. output/pipeline_info/MetaAnalysis_report.html.
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from work folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - .command.sh: script which was submitted
  - .command.log: log file for seeing the analysis outputs/errors.
  - .command.err: file which lists the errors, if any.

**NB!** Pipeline does not overwrite the files in the output folder. Therefore, **before re-running the pipeline, delete the output files in the output folder!**

##### Output

Pipeline produces many parquet (`.parquet`) files to the output directory which contain the summary statistics from genome-wide eQTL analysis. The output should take ~6TB of disk space.

Those numerous `.parquet` files can now further processed:

- To extract all the effects which reach to certain significance level.
- To extract the genome-wide association profile(s) for certain genes of interest.

These tasks can be achieved by separate pipeline [ExtractMetaAnalysisResults](https://github.com/eQTLGen/ExtractMetaAnalysisResults)

#### Acknowledgements

This pipeline utilizes HASE (https://github.com/roshchupkin/hase) and some of its helper scripts originally developed by:

Gennady V. Roscupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands)

Hieab H. Adams (Department of Epidemiology, Erasmus MC, Rotterdam, Netherlands).

##### Changes from the original HASE repo

Robert Warmerdam (Department of Genetics, University Medical Center Groningen, University of Groningen, Groningen, Netherlands) modified the original HASE and fixed some bugs.

Urmo Võsa (Institute of Genomics, University of Tartu, Tartu, Estonia) incorporated it into Nextflow pipeline and customized/supplanted some parts of the code.

Changes:

- Fixed bug causing an exception when more than 1000 individuals were used.
- Resolved bug causing the --intercept option having no effect.
- Made version numbers of pip packages explicit.
- Added commentary to code in places.
- Setting /hdgwas/data.py were changed self.chunk_size=10000 --> self.chunk_size=20000.
- Use PyCryptodome for encoding datasets.
- Implemented inverse variance meta-analysis ion addition to original combined mega-analysis.
- Implemented an enhancement which allows missing values in meta-analysis (e.g. when some gene or variant fails QC in some specific dataset then effects from remaining datasets are used).
- Partial derivative matrices are read into memory by chunk.
- For eQTLGen pipelines: removed folders with unit tests and test data, in order to keep the tool lightweight.

##### Citation

Original method paper for HASE framework:

Roshchupkin, G. V. et al. HASE: Framework for efficient high-dimensional association analyses. Sci. Rep. 6, 36076; doi: 10.1038/srep36076 (2016)
