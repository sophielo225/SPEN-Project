# This file contains codes to make tables

# Load all required library
library(tidyverse)
library(knitr)

# Make tables for tsv files acquired from GSEA website

# Helper function that saves tsv files as Markdown texts
save_tsv_as_md <- function(file_path, out_dir = "tables") {
    
    if (!dir.exists(out_dir)) dir.create(out_dir)
    
    short_name <- tools::file_path_sans_ext(basename(file_path))
    
    lines <- readLines(file_path)
    start <- grep("Gene/Gene Set Overlap Matrix", lines)
    matrix_lines <- lines[(start + 2):length(lines)]
    matrix_lines <- matrix_lines[matrix_lines != ""]
    
    df <- read.delim(text = matrix_lines, sep = "\t", header = TRUE)
    
    # Modify the data frame
    df <- df %>%
        dplyr::select(-`Gene.Description`) %>%
        dplyr::rename(Entrez_ID = `Entrez.Gene.Id`) %>%
        dplyr::rename(Gene_symbol = `Gene.Symbol`) %>%
        pivot_longer(cols = -c(Entrez_ID, Gene_symbol), names_to = "pathway_type", values_to = "pathway") %>%
        filter(pathway != "" & !is.na(pathway)) %>%
        group_by(Entrez_ID, Gene_symbol) %>%
        dplyr::summarize(Pathway = paste(sort(unique(pathway)), collapse = ", "), .groups = "drop")
    
    kable(df, format = "simple") %>%
        writeLines(paste0(out_dir, "/", short_name, "_pathway_results.md"))
}

# Apply function to all tsv files
tsv_files <- list.files("pathway_results", pattern = "\\.tsv$", full.names = TRUE)
lapply(tsv_files, save_tsv_as_md)