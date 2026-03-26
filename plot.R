# This file contains codes to make plots or tables

# Load all required library
library(EnhancedVolcano)
library(tidyverse)

# Graph volcano plots for GSE43795 and E-MEXP-1914
# Load RDS tibbles
gse43795_results <- readRDS("gse43795_results_tibble.rds")
emexp_1914_results <- readRDS("emexp_1914_results_tibble.rds")

top_genes <- rownames(gse43795_results)[order(gse43795_results$adj.P.Val)][1:10] # Only label top genes

pdf("Volcano_plot_GSE43795.pdf", width = 10, height = 8)
EnhancedVolcano(gse43795_results,
                lab = rownames(gse43795_results),
                x = 'logFC',
                y = 'adj.P.Val',
                selectLab = top_genes,
                
                pCutoff = 0.05,
                FCcutoff = 1,
                
                pointSize = 1.5,
                labSize = 3,
                colAlpha = 0.5,
                col = c("grey70", "#56B4E9", "#56B4E9", "#D55E00"),
                # ylim = c(0, 12),
                
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 100,
                
                title = 'Volcano Plot for GSE43795',
                subtitle = 'Differential Expression',
                caption = 'log2FC cutoff = 1, FDR < 0.05'
)
dev.off()

# Make tables for tsv files acquired from GSEA website
data <- read_tsv("GSEA_result_10_genes.tsv")
