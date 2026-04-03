# This file contains codes to make plots or tables

# Load all required library
library(EnhancedVolcano)
library(tidyverse)
library(knitr)

#####################################
# Graph volcano plots for GSE43795 and E-MEXP-1914

# Load RDS tibbles
gse43795_results <- readRDS("gse43795_results_tibble.rds")
emexp_1914_results <- readRDS("emexp_1914_results_tibble.rds")

# Create folder for plots
if (!dir.exists("figures")) {
    dir.create("figures")
}

# Plot GSE43795
gse43795_top_genes <- gse43795_results %>%  # Only label top 10 genes by p-value
    arrange(adj.P.Val) %>%
    dplyr::slice(1:10) %>%
    pull(hgnc_symbol)

pdf("figures/Volcano_plot_GSE43795.pdf", width = 10, height = 8)
EnhancedVolcano(gse43795_results,
                lab = gse43795_results$hgnc_symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                selectLab = gse43795_top_genes,
                
                pCutoff = 0.05,
                FCcutoff = 1,
                
                pointSize = 1.5,
                labSize = 3,
                colAlpha = 0.2,
                col = c("grey70", "#56B4E9", "#56B4E9", "#D55E00"),
                
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 100,
                
                title = 'Volcano Plot for GSE43795',
                subtitle = 'Differential Expression',
                caption = 'log2FC cutoff = 1, FDR < 0.05'
)
dev.off()

# Plot E-MEXP-1914
emexp_1914_top_genes <- emexp_1914_results %>%  # Only label top 10 genes by p-value
    arrange(adj.P.Val) %>%
    dplyr::slice(1:10) %>%
    pull(hgnc_symbol)

pdf("figures/Volcano_plot_EMEXP_1914.pdf", width = 10, height = 8)
EnhancedVolcano(emexp_1914_results,
                lab = emexp_1914_results$hgnc_symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                selectLab = emexp_1914_top_genes,
                
                pCutoff = 0.05,
                FCcutoff = 1,
                
                pointSize = 1.5,
                labSize = 3,
                colAlpha = 0.2,
                col = c("grey70", "#56B4E9", "#56B4E9", "#D55E00"),
                
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 100,
                
                title = 'Volcano Plot for E-MEXP-1914',
                subtitle = 'Differential Expression',
                caption = 'log2FC cutoff = 1, FDR < 0.05'
)
dev.off()

#####################################
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