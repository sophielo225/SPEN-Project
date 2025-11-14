# This file analyzes two mRNA sequencing data sets(GSE43795, E-MEXP-1914)
# GSE43795 is from GEO, and E-MEXP-1914 is from ArrayExpress

# Install needed packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("ArrayExpress")

# Load libraries
library(GEOquery)
library(tidyverse)
library(ArrayExpress)

# Get expression data, metadata, and feature data
# GSE43795
gse43795 = getGEO("GSE43795", GSEMatrix = TRUE)

gse43795_expression_data = exprs(gse43795[[1]]) %>% 
    as_tibble(rownames = "ID_REF")

gse43795_metadata = pData(gse43795[[1]]) %>%
    as_tibble(rownames = "Sample_ID")

gse43795_feature_data = fData(gse43795[[1]]) %>%
    as_tibble(rownames = "ID_REF")

# E-MEXP-1914
ae <- ArrayExpress("E-MEXP-1914", fix = TRUE)
emexp1914_expression_data = exprs(ae)
emexp1914_metadata = pData(ae)
emexp1914_feature_data = fData(ae)
