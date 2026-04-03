# This file contains codes to make plots

# Load all required library
library(EnhancedVolcano)
library(tidyverse)

# Graph volcano plots for GSE43795 and E-MEXP-1914

# Load RDS tibbles
gse43795_results <- readRDS("gse43795_results_tibble.rds")
emexp_1914_results <- readRDS("emexp_1914_results_tibble.rds")

# Create folder for plots
if (!dir.exists("figures")) {
    dir.create("figures")
}

# Plot GSE43795
# Rank genes by adjusted p-values and log fold change
gse43795_top_genes_df <- gse43795_results %>% 
    arrange(adj.P.Val) %>%
    mutate(rank_pvalue = row_number()) %>%  
    arrange(desc(abs(logFC))) %>%
    mutate(rank_logFC = row_number()) %>%
    mutate(rank_avg = (rank_pvalue + rank_logFC) / 2) %>%
    arrange(rank_avg)

gse43795_top_genes <- gse43795_top_genes_df %>%
    dplyr::slice(1:15) %>%  # Only label top 15 genes
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
# Rank genes by adjusted p-values and log fold change
emexp_1914_top_genes_df <- emexp_1914_results %>% 
    arrange(adj.P.Val) %>%
    mutate(rank_pvalue = row_number()) %>%  
    arrange(desc(abs(logFC))) %>%
    mutate(rank_logFC = row_number()) %>%
    mutate(rank_avg = (rank_pvalue + rank_logFC) / 2) %>%
    arrange(rank_avg)

emexp_1914_top_genes <- emexp_1914_top_genes_df %>%
    dplyr::slice(1:15) %>%  # Only label top 15 genes
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
