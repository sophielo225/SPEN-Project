# Install all required library
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("GEOquery","oligo"), ask = FALSE)

if (!requireNamespace("R.utils", quietly = TRUE))
    install.packages("R.utils")

library(GEOquery)
library(oligo)
library(R.utils)
library(pd.mirna.4.0)
library(affyio)
library(tidyverse)
library(broom)

gse_id = "GSE140719"
cel_dir = file.path(gse_id, "CEL")

# Unzip all .cel.gz files
cel_archives = list.files(cel_dir, pattern = "\\.cel\\.gz$", full.names = TRUE)
for (f in cel_archives) {
    message("Unzipping: ", basename(f))
    R.utils::gunzip(f, overwrite = TRUE, remove = TRUE)  # remove the .gz after extraction
}

# Check results. It should show .cel files 
list.files(cel_dir)


cel_files <- list.celfiles(cel_dir, full.names = TRUE)
stopifnot(length(cel_files) > 0)

# Read raw CELs without annotation
raw <- read.celfiles(cel_files, pkgname = "pd.mirna.4.0")

# Background correction + quantile normalization + summarization
eset <- rma(raw)                  # ExpressionSet
expr_mat <- Biobase::exprs(eset)  # numeric matrix (features x samples)

# Save results
out_file <- file.path(gse_id, paste0(gse_id, "_RMA_expression.tsv"))
write.table(expr_mat, out_file, sep = "\t", quote = FALSE, col.names = NA)

message("RMA matrix written to: ", out_file)

# Now we have the normalized expression dataset
# Read it into tibble and modify column names
gse140719_expression_data = read_tsv("GSE140719/GSE140719_RMA_expression.tsv")
gse140719_expression_data = rename(gse140719_expression_data, ID_REF = `...1`)
colnames(gse140719_expression_data) = sub("_.*\\.cel$", "", colnames(gse140719_expression_data))

# Get metadata and feature data
gse140719 = getGEO("GSE140719", GSEMatrix = TRUE)

# gse140719_expression_data = exprs(gse140719[[1]]) %>%
#     as_tibble(rownames = "ID_REF")

gse140719_metadata = pData(gse140719[[1]]) %>%
    as_tibble(rownames = "Sample_ID")

gse140719_feature_data = fData(gse140719[[1]]) %>%
    as_tibble(rownames = "ID_REF")

# Filter only miRNA entries and map MIMAT ID, transcript ID to expression data
gse140719_MIMAT_ID_REF = rename(gse140719_feature_data, Sequence_Type = 'Sequence Type') %>%
    rename(Transcript_ID = `Transcript ID(Array Design)`) %>%
    filter(Sequence_Type == "miRNA") %>%
    select(ID_REF, Accession, Transcript_ID)

gse140719_expression_data = left_join(gse140719_expression_data, gse140719_MIMAT_ID_REF, by = "ID_REF") %>%
    rename(MIMAT_ID = Accession) %>%
    select(MIMAT_ID, Transcript_ID, GSM4182441:GSM4182451) %>%
    filter(grepl("^MIMAT", MIMAT_ID))

# Prepare data to do ANOVA
gse140719_samples = select(gse140719_metadata, Sample_ID, title) %>%
    mutate(Group = case_when(
        grepl("Normal", title) ~ "normal",
        grepl("LT", title) ~ "localized",
        grepl("MT", title) ~ "metastatic")) %>%
    select(Sample_ID, Group)

# Map groups into expression data
gse140719_three_groups = pivot_longer(gse140719_expression_data, cols = starts_with("GSM"),
                                      names_to = "Sample_ID",values_to = "Expression") %>%
    left_join(gse140719_samples, by = "Sample_ID")

# Perform ANOVA
anova_results <- gse140719_three_groups %>%
    group_by(MIMAT_ID, Transcript_ID) %>%  # group by each miRNA
    do(tidy(aov(Expression ~ Group, data = .))) %>%  # run ANOVA within each group
    filter(term == "Group") %>%  # only keep the main Group effect
    select(MIMAT_ID, Transcript_ID, p.value, statistic, df, sumsq, meansq)  # optional: keep key columns

# Perform Benjamini-Hochberg correction
anova_results <- anova_results %>%
    mutate(p_adj = p.adjust(p.value, method = "BH"))

anova_results %>%
    filter(p_adj < 0.05) %>%
    arrange(p_adj) -> anova_significant 

# Filter the data set to get significant IDs and calculate mean expression per group
anova_significant %>% pull(MIMAT_ID) -> significant_MIMAT_ID
gse140719_three_groups %>%
    filter(MIMAT_ID %in% significant_MIMAT_ID) -> significant_expr_data

significant_expr_data %>%
    group_by(MIMAT_ID, Transcript_ID, Group) %>%
    summarise(mean_expr = mean(Expression)) %>%
    pivot_wider(names_from = Group, values_from = mean_expr) -> mean_expr_data

# Do fold change calculation
mean_expr_data %>%
    mutate(FC_localized = localized - normal) %>%
    mutate(FC_metastatic = metastatic - normal) %>%
    mutate(FC_localized_vs_metastatic = localized - metastatic) -> fold_change_data

# Apply fold change cutoff (>2 or <-2)
# Since the data has been log2 transformed, so the fold change cutoff would be 1 
fold_change_data %>%
    filter(abs(FC_localized_vs_metastatic) > 1) %>%
    pull(MIMAT_ID) -> significant_MIMAT_ID
    
anova_results %>%
    filter(MIMAT_ID %in% significant_MIMAT_ID) %>%
    arrange(p_adj)-> eleven_IDs


# Visualize the miRNA 
gse140719_three_groups %>%
    filter(Transcript_ID == "hsa-miR-4448") %>%
    ggplot(aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot() +
    theme_bw()
    labs(title = "Expression of hsa-miR-4448 by Group")

    

# Perform Mann-Whitney U test for localized vs metastatic
# Prepare data to do Mann-Whitney U test
gse140719_samples = gse140719_metadata %>%
    filter(grepl("LT-localized tumor", title) | grepl("MT-Metastatic tumor", title)) %>%
    mutate(Group = if_else(grepl("LT-localized tumor", title), "localized", "metastatic")) %>%
    select(Sample_ID, Group)

# Map groups into expression data
gse140719_two_groups = pivot_longer(gse140719_expression_data, cols = starts_with("GSM"),
                                      names_to = "Sample_ID",values_to = "Expression") %>%
    left_join(gse140719_samples, by = "Sample_ID") %>%
    filter(Group == "localized" | Group == "metastatic")

# Perform Mann-Whitney U test
gse140719_result = gse140719_two_groups %>%
    group_by(MIMAT_ID, Transcript_ID) %>%
    summarize(p_value = wilcox.test(Expression ~ Group)$p.value) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))

# Get significant IDs
gse140719_result %>%
    filter(p_adj < 0.05) %>%
    arrange(p_adj) -> mann_whitney_significant 

# Filter the data set to get significant IDs and calculate mean expression per group
mann_whitney_significant %>% pull(MIMAT_ID) -> significant_MIMAT_ID
gse140719_two_groups %>%
    filter(MIMAT_ID %in% significant_MIMAT_ID) -> significant_expr_data

significant_expr_data %>%
    group_by(MIMAT_ID, Transcript_ID, Group) %>%
    summarise(mean_expr = mean(Expression)) %>%
    pivot_wider(names_from = Group, values_from = mean_expr) -> mean_expr_data

# Do fold change calculation
mean_expr_data %>%
    mutate(fold_change = localized - metastatic) -> fold_change_data

# Apply fold change cutoff (>2 or <-2)
fold_change_data %>%
    filter(abs(fold_change) > 1) %>%
    pull(Transcript_ID) %>%
    print()
