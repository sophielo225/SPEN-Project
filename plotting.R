# This file contains codes to make plots or tables

# Load all required library
library(EnhancedVolcano)
library(tidyverse)

# Graph volcano plot

# Make tables for tsv files acquired from GSEA website
data <- read_tsv("GSEA_result_10_genes.tsv")
