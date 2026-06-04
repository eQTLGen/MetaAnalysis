#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(readxl)
library(data.table)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  # After writing this file, use the following in R

  explained_variance <- fread("explained_variance_2024-06-26.txt")
  cohort_mapping <- explained_variance %>% distinct(cohort, cohort_new_name)
  covariates <- read_excel("covariates_long_2024-06-23.xlsx", sheet="covariates_long") %>%
    rename(cohort_new_name = "cohort") %>%
    inner_join(cohort_mapping, unmatched = "error") %>%
    left_join(explained_variance, unmatched = "error")

  permuted <- fread("mastertable_permuted_2024-06-26_extended_fixedSampleNames_fixedGeneIds.txt") %>%
    inner_join(covariates, unmatched = "error")
  empirical <- fread("mastertable_empirical_2024-05-17_extended.txt") %>%
    inner_join(covariates, unmatched = "error")

  covariates_extended <- covariates %>%
    mutate(
      analysis2_4GenPC20ExpPC = case_when(analysis1_4GenPCNoExpPC | PC %in% paste0("PC", 1:20) ~ T, TRUE ~ F),
      analysis6_4GenPC50ExpPC = case_when(analysis1_4GenPCNoExpPC | PC %in% paste0("PC", 1:50) ~ T, TRUE ~ F),
      analysis3_4GenPC25pExpPC = case_when(analysis1_4GenPCNoExpPC | cumulative_explained_variance < 0.25 | lag(cumulative_explained_variance) < 0.25 ~ T, TRUE ~ F),
      analysis4_4GenPC50pExpPC = case_when(analysis1_4GenPCNoExpPC | cumulative_explained_variance < 0.50 | lag(cumulative_explained_variance) < 0.50 ~ T, TRUE ~ F),
      analysis5_4GenPC60pExpPC = case_when(analysis1_4GenPCNoExpPC | cumulative_explained_variance < 0.60 | lag(cumulative_explained_variance) < 0.60 ~ T, TRUE ~ F),
    ) %>%
    select(cohort_new_name, cohort, covariate, PC, comment, explained_variance, cumulative_explained_variance, analysis1_4GenPCNoExpPC, analysis2_4GenPC20ExpPC, analysis3_4GenPC25pExpPC, analysis4_4GenPC50pExpPC, analysis5_4GenPC60pExpPC, analysis6_4GenPC50ExpPC)

  fwrite(covariates_extended, "~/eQTLGen/freeze3/eqtl_mapping/input/covariates_long_extended_2024-08-08.txt", col.names=T, row.names=F, quote=F, sep="\t", na="NA")

  fwrite(table %>% filter(analysis2_4GenPCNoExpPC) %>%
  group_by(cohort) %>%
  summarise(covariate_list=list(covariate)) %>%
  unnest_wider(covariate_list)


  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}