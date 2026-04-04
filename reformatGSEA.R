# This file contains codes to make tables

# Load all required library
library(tidyverse)
library(knitr)

# Make tables for tsv files acquired from GSEA website

# Helper function that saves tsv files as Markdown texts
save_tsv_as_md <- function(file_path, out_dir = "tables") {
    
    if (!dir.exists(out_dir)) dir.create(out_dir)
    
    short_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Read lines
    lines <- readLines(file_path)
    result_lines <- lines[10:20]
    result_df <- read_tsv(paste(result_lines, collapse = "\n"))
    matrix_df <- read_tsv(file_path, skip = 24)
    
    # Convert overlap matrix to long format
    matrix_long <- matrix_df %>%
        pivot_longer(cols = -(c(`Entrez Gene Id`, `Gene Symbol`, `Gene Description`)),
            names_to = "GeneSet", values_to = "Present") %>%
        filter(!is.na(Present) & Present != "")
    
    # Collapse gene symbols per gene set
    gene_lists <- matrix_long %>%
        group_by(GeneSet) %>%
        dplyr::summarize(Gene_Symbols = paste(`Gene Symbol`, collapse = ", "), .groups = "drop")
    
    # Join with overlap summary table
    final_table <- result_df %>%
        dplyr::select(`Gene Set Name`, `p-value`, `FDR q-value`) %>%
        left_join(gene_lists, by = c("Gene Set Name" = "GeneSet")) %>%
        dplyr::rename(Gene_Set_Name = `Gene Set Name`) %>%
        dplyr::rename(FDR_q_value = `FDR q-value`) %>%
        dplyr::rename(p_value = `p-value`)
    
    # Write the table into markdown text
    kable(final_table, format = "simple", digits = 30) %>%
        writeLines(paste0(out_dir, "/", short_name, "_pathway_results.md"))
}

# Apply function to all tsv files
tsv_files <- list.files("pathway_results", pattern = "\\.tsv$", full.names = TRUE)
lapply(tsv_files, save_tsv_as_md)
