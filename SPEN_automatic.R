# Install all required library

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)

# Download GEO Series
gse43796 = getGEO("GSE43796", GSEMatrix = TRUE)
gse140719 = getGEO("GSE140719", GSEMatrix = TRUE)


# Get GSE43796 expression data, metadata, and feature data (probe annotation)
gse43796_expression_data = exprs(gse43796[[1]]) %>% 
    as_tibble(rownames = "ID_REF")