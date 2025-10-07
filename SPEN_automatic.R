# Install all required library

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)

# Download GEO Series
gse43796 = getGEO("GSE43796", GSEMatrix = TRUE)
gse140719 = getGEO("GSE140719", GSEMatrix = TRUE)

# GSE43796
# Get GSE43796 expression data, metadata, and feature data (probe annotation)
gse43796_expression_data = exprs(gse43796[[1]]) %>% 
    as_tibble(rownames = "ID_REF")

gse43796_metadata <- pData(gse43796[[1]]) %>%
    as_tibble(rownames = "Sample_ID")

gse43796_feature_data <- fData(gse43796[[1]]) %>%
    as_tibble(rownames = "ID_REF")

# Get MIMAT ID from feature data and map MIMAT ID to expression data
gse43796_MIMAT_ID_REF = mutate(gse43796_feature_data, MIMAT_ID = str_extract(ACCESSION_STRING, "MIMAT\\d+")) %>%
    select(ID_REF, MIMAT_ID)
gse43796_expression_data = left_join(gse43796_expression_data, gse43796_MIMAT_ID_REF, by = "ID_REF") %>%
    select(MIMAT_ID, GSM1071387:GSM1071417)

# Since GSE43796 has duplicate expression level for some miRNA ID, 
# we need to pivot the tibble and calculate the average level if the miRNA ID has multiple expression levels
gse43796_expression_data = pivot_longer(gse43796_expression_data, GSM1071387:GSM1071417, names_to = "Sample_name", values_to = "Expression_level") %>%
    group_by(MIMAT_ID, Sample_name) %>%
    summarize(mean_expression = mean(Expression_level)) %>%
    pivot_wider(names_from = "Sample_name", values_from = "mean_expression")

# Get normal and tumor samples from metadata
gse43796_sample_group = gse43796_metadata %>% 
    filter(source_name_ch1 == "solid-pseudopapillary neoplasm" | source_name_ch1 == "non-neoplastic pancreas") %>%
    mutate(Group = if_else(source_name_ch1 == "solid-pseudopapillary neoplasm", "tumor", "normal")) %>%
    select(Sample_ID, Group)

# Map sample group to expression data
gse43796_clean = pivot_longer(gse43796_expression_data, cols = starts_with("GSM"),
    names_to = "Sample_ID",values_to = "Expression") %>%
    left_join(gse43796_sample_group, by = "Sample_ID") %>%
    filter(Group == "tumor" | Group == "normal")


# GSE140719
# Get GSE140719 expression data, metadata, and feature data (probe annotation)
gse140719_expression_data = exprs(gse140719[[1]]) %>% 
    as_tibble(rownames = "ID_REF")

gse140719_metadata <- pData(gse140719[[1]]) %>%
    as_tibble(rownames = "Sample_ID")

gse140719_feature_data <- fData(gse140719[[1]]) %>%
    as_tibble(rownames = "ID_REF")

# Filter only miRNA entries and map MIMAT ID to expression data
gse140719_MIMAT_ID_REF = rename(gse140719_feature_data, Sequence_Type = 'Sequence Type') %>%
    filter(Sequence_Type == "miRNA") %>%
    select(ID_REF, Accession)
gse140719_expression_data = left_join(gse140719_expression_data, gse140719_MIMAT_ID_REF, by = "ID_REF") %>%
    filter(grepl("^MIMAT", ID_REF)) %>%
    rename(MIMAT_ID = Accession) %>%
    select(MIMAT_ID, GSM4182441:GSM4182451)

# Get normal and tumor samples from metadata
gse140719_sample_group = gse140719_metadata %>% 
    filter(grepl("LT-localized tumor", title) | grepl("Normal", title)) %>%
    mutate(Group = if_else(grepl("LT-localized tumor", title), "tumor", "normal")) %>%
    select(Sample_ID, Group)

# Map sample group to expression data
gse140719_clean = pivot_longer(gse140719_expression_data, cols = starts_with("GSM"),
                              names_to = "Sample_ID",values_to = "Expression") %>%
    left_join(gse140719_sample_group, by = "Sample_ID") %>%
    filter(Group == "tumor" | Group == "normal")

# Select mutual MIMAT ID from both data sets
gse43796_MIMAT = pull(gse140719_expression_data, MIMAT_ID)
gse140719_MIMAT = pull(gse43796_expression_data, MIMAT_ID)
mutual_MIMAT = intersect(gse43796_MIMAT, gse140719_MIMAT)
gse43796_clean = filter(gse43796_clean, MIMAT_ID %in% mutual_MIMAT)
gse140719_clean = filter(gse140719_clean, MIMAT_ID %in% mutual_MIMAT)


# Perform Mann-Whitney Test and find q-value
gse43796_result = gse43796_clean %>%
    group_by(MIMAT_ID) %>%
    summarize(p_value = wilcox.test(Expression ~ Group)$p.value) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))

gse140719_result = gse140719_clean %>%
    group_by(MIMAT_ID) %>%
    summarize(p_value = wilcox.test(Expression ~ Group)$p.value) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))

# Find statistically significant MIMAT IDs that are in both data sets
# Set the threshold into 0.05
gse43796_result %>%
    filter(q_value < 0.05) %>%
    pull(MIMAT_ID) -> gse43796_MIMAT_sig
print(gse43796_MIMAT_sig)

gse140719_result %>%
    filter(q_value < 0.05) %>%
    pull(MIMAT_ID) -> gse140719_MIMAT_sig
print(gse140719_MIMAT_sig)
