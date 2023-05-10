library(data.table)
library(stringr)
library(fs)
library(tidyr)
library(gtools)

setDTthreads(8)

#########################
# arguments and settings#
#########################

outputpath <- "[path to your output]"

dataset_path <- c("/[path to cohort 1]/[name of cohort 1]_IntermediateFilesEncoded_to_upload/",
"/[path to cohort 2]/[name of cohort 2]_IntermediateFilesEncoded_to_upload/",
"/[path to cohort 3]/[name of cohort 3]_IntermediateFilesEncoded_to_upload/")

# You can adjust the SNP QC thresholds here
MAF_thresh <- 0.01
HWE_thresh <- 1e-6
MACH_R2_thresh <- 0.4
CALL_thresh <- 1

MAC_thresh <- 10

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

# Make SNP inclusion files
fileConn<-file(paste0(outputpath, "/covariate_indices.txt"), "w")

sink(paste0(outputpath, "/filter_logs.log"))
cat("Dataset\tN\tNr. SNPs\tFiltered SNPs\tTyped SNPs\tSNPs MAC>10; not filtered\tSNPs AC>5; not filtered\tGenes original\tGenes filtered\tGenetically regulated PCs\n")

for (dataset in dataset_path){

dataset_name <- str_replace(as.character(fs::dir_ls(path = paste0(dataset, "/SumStats/genotypes/SNPQC/"), glob = "*_SNPQC.txt.gz")), ".*/SNPQC/", "")
dataset_name <- str_replace(dataset_name, "_SNPQC.txt.gz", "")

message(paste0("Processing ", dataset_name))
message("Make SNP inclusion files")

snpqc_filename <- fs::dir_ls(path = paste0(dataset, "/SumStats/genotypes/SNPQC/"), glob = "*_SNPQC.txt.gz")
snpqc <- fread(snpqc_filename)

# fix initial filter MAF>0.01
snpqc <- snpqc[MAF > 0.01]

# filters
orig_snps <- nrow(snpqc)
snpqc <- snpqc[MAF > MAF_thresh & HWE > HWE_thresh & MACH_R2 > MACH_R2_thresh & CALL == CALL_thresh]

# Count variants which have at least 5 allele counts for rarer genotype
hom1_gen_counts <- str_replace(snpqc$Genotype_Counts, ", .*", "")
snpqc$hom1_gen_counts <- as.numeric(str_replace(hom1_gen_counts, ".*: ", ""))

hom2_gen_counts <- str_replace(snpqc$Genotype_Counts, ".*: ", "")
snpqc$hom2_gen_counts <- as.numeric(str_replace(hom2_gen_counts, ",", ""))

het_gen_counts <- sub('(,.*?),.*', '\\1', snpqc$Genotype_Counts)
snpqc$het_gen_counts <- as.numeric(str_replace(het_gen_counts, ".* ", ""))

# Calculate MAC for each SNP
mac <- snpqc[, .(MAC = 2 * min(c(hom1_gen_counts, hom2_gen_counts)) + het_gen_counts), by = .(ID)]

# Calculate N
N_split <- unlist(str_split(snpqc$Genotype_Counts[1], ", "))
N_split <- str_replace(N_split, ".*: ", "")
N_split <- as.numeric(str_replace(N_split, ",", ""))
N_split <- sum(N_split)

snpqc <- snpqc[ID %in% mac[MAC > MAC_thresh]$ID]

fwrite(data.table(ID = snpqc$ID), paste0(outputpath, "/", dataset_name, "_SnpsToInclude.txt"), sep = "\t", row.names = FALSE)

geneqc_processed <- fread(paste0(dataset, "/SumStats/expression/processed_gene_summary.txt"))
genes_to_include <- geneqc_processed[geneqc_processed$nr_unique_values > gene_uniq_thresh * max(geneqc_processed$nr_values)]$gene
fwrite(data.frame(ID = genes_to_include), paste0(outputpath, "/", dataset_name, "_GenesToInclude.txt"), sep = "\t", row.names = FALSE)

message("Check genetically regulated expression PCs")
# Report PCs under genetic regulation, make covariate indices file
genetic_pcs <- fread(paste0(dataset, "/", dataset_name, "_GenRegPcs.txt.gz"))
genetic_pcs <- genetic_pcs[`P-value` < GenRegPC_thresh]
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
nrow(snpqc[snpqc$TYPED == "1", ]),
nrow(mac[MAC > 10]),
nrow(snpqc[hom1_gen_counts > 5 & hom2_gen_counts > 5]),
nrow(geneqc_processed),
length(genes_to_include),
paste(genetic_pcs, collapse = ", "),
"\n", sep = "\t"))

}

message("Finalising...")
sink()