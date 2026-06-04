#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)
library(optparse)

# Declare constants
option_list <- list(
  make_option(c("--mastertable"), type="character", help="Mastertable to check data for"),
)

parser <- OptionParser(usage="%prog [options]", description="Validate master table", option_list=option_list)



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
  mastertable <- fread(args$mastertable)

  # Ho

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}