#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)
library(optparse)

# Declare constants
option_list <- list(
  make_option(c("--file_name_mapping"), type="character", help="File with cohort names and cohort locations to use"),
  make_option(c("--encoded"), type="character", default="1", help="Is there individual level data or encoded data. Assumed to be constant over all datasets."),
  make_option(c("--permuted"), action="store_true", default=FALSE, help="Generate master table for permuted data."),
  make_option(c("--inclusion_list_dir"), type="character", help="Directory where inclusion lists are stored."),
  make_option(c("--out"), type="character", help="Output file")
)

parser <- OptionParser(usage="%prog [options]", description="Generate Master Table from fileNameMapping", option_list=option_list)



# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  args <- parser$parse(argv)

  # Process input
  file_name_mapping <- fread(args$file_name_mapping)

  # Perform method
  master_table <- file_name_mapping %>%
    mutate(encoded = 1,
           genotype = file.path(dataset_path, "empirical", "EncodedGenotypeData"),
           expression = file.path(dataset_path, "empirical", "EncodedPhenotypeData"),
           partial_derivatives = file.path(dataset_path, "empirical", "pd_shared"),
           snp_inclusion = file.path(inclusion_list_dir, paste0(dataset_name_new, "_SnpsToInclude.txt")),
           gene_inclusion = file.path(inclusion_list_dir, paste0(dataset_name_new, "_GenesToInclude.txt"))) %>%
    select(c("cohort" = "dataset_name", "cohort_new_name" = "dataset_name_new"), "genotype", "expression", "partial_derivatives", "encoded", "snp_inclusion", "gene_inclusion")

  fwrite(master_table, args$out, sep="\t", row.names=F, col.names=T, quote=F)
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}