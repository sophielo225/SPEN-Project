# Load all the required libraries
library(GEOquery)
library(tidyverse)
library(limma)

# GSE43796
# Get metadata and feature data
gse43796 = getGEO("GSE43796", GSEMatrix = TRUE)

gse43796_metadata = pData(gse43796[[1]]) %>%
    as_tibble(rownames = "Sample_ID") %>%
    unite(Sample_ID, Sample_ID, description, sep="_") %>%
    filter(source_name_ch1 == "solid-pseudopapillary neoplasm" | source_name_ch1 == "non-neoplastic pancreas")

gse43796_feature_data = fData(gse43796[[1]]) %>%
    as_tibble(rownames = "ID_REF")

# Download all supplementary files for the series (includes the .tar raw file)
getGEOSuppFiles("GSE43796")

# Set the file path
tarfile <- "GSE43796/GSE43796_RAW.tar"

# Unpack the tar file into a subfolder
untar(tarfile, exdir = "GSE43796/raw_data")

# Gunzip them if they are compressed
txt_files <- list.files("GSE43796/raw_data", pattern = "\\.gz$", full.names = TRUE)
sapply(txt_files, R.utils::gunzip, overwrite=TRUE)

raw_files <- list.files("GSE43796/raw_data", pattern = "\\.txt$", full.names = TRUE)

# Read Agilent one-color arrays
RG <- read.maimages(raw_files, source="agilent", green.only=TRUE)

# Do normalization
# 'normexp' is widely used and similar in concept to RMA background correction
RG_bc <- backgroundCorrect(RG, method="normexp")  

# Normalize between arrays
MA <- normalizeBetweenArrays(RG_bc, method="quantile")

# Remove directory names from the column names in the expression matrix
colnames(MA$E) <- gsub("GSE43796/raw_data/", "", colnames(MA$E))

# Only keep the samples we care about
MA <- MA[,pull(gse43796_metadata, Sample_ID)]

# Only keep the regular probes (remove control probes)
MA <- MA[MA$genes$ControlType==0, ]

# Map MIMAT IDs to each probe names in MA$genes
gse43796_MIMAT_ID_REF <- mutate(gse43796_feature_data, MIMAT_ID = str_extract(ACCESSION_STRING, "MIMAT\\d+")) %>%
    dplyr::select(ID_REF, MIMAT_ID) %>%
    dplyr::rename(ProbeName = ID_REF)

MA$genes <- MA$genes %>%
    left_join(gse43796_MIMAT_ID_REF, by = "ProbeName")

# Get the mean expression for each MIMAT ID
MA_collapse <- avereps(MA, ID=MA$genes$MIMAT_ID)

# Create the sample group factor
samples <- colnames(MA_collapse$E)
group <- ifelse(grepl("_S-", samples), "tumor", "normal")
group <- factor(group, levels=c("normal","tumor"))

# Build the design matrix
design <- model.matrix(~group)

# Fit the linear model and find DE (differentially expressed) probes
fit <- lmFit(MA_collapse, design)
fit <- eBayes(fit, trend=TRUE)
summary(decideTests(fit))
top <- topTable(fit, coef="grouptumor", number=Inf, adjust.method="BH")

# Convert the result matrix into tibble
gse43796_result <- as_tibble(top) %>%
    dplyr::select(MIMAT_ID, SystematicName, logFC, adj.P.Val)

# Get the significant IDs (adjusted p-value < 0.05 and fold change > 2)
gse43796_significant = gse43796_result %>%
    filter(abs(logFC) > 1 & adj.P.Val < 0.05)


#######################################
# GSE140719 Analysis
# Load all libraries
library(oligo)
library(R.utils)
library(pd.mirna.4.0)
library(affyio)

# Get metadata and feature data
gse140719 = getGEO("GSE140719", GSEMatrix = TRUE)

gse140719_metadata = pData(gse140719[[1]]) %>%
    as_tibble(rownames = "Sample_ID") %>%
    filter(grepl("localized", title) | grepl("Normal", title)) %>%
    dplyr::select(Sample_ID, title)

gse140719_feature_data = fData(gse140719[[1]]) %>%
    as_tibble(rownames = "ID_REF") %>%
    dplyr::select(ID_REF, Accession)

# Get CEL files
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

# Tidy column names of the expression matrix
colnames(expr_mat) <- gsub("_.*", "", colnames(expr_mat))

# Only get the samples we are interested in the expression matrix
expr_mat <- expr_mat[,pull(gse140719_metadata, Sample_ID)]

# Map accession IDs to expression matrix
probe_tbl <- tibble(ID_REF = rownames(expr_mat))
probe_tbl <- left_join(probe_tbl, gse140719_feature_data, by = "ID_REF")
new_rownames = probe_tbl$Accession
rownames(expr_mat) <- new_rownames
expr_mat_mimat <- expr_mat[grepl("MIMAT", rownames(expr_mat)), ]

# Create the sample group factor
gse140719_group <- ifelse(grepl("Normal", gse140719_metadata$title),
                "normal", "tumor")
gse140719_group <- factor(gse140719_group, levels = c("normal", "tumor"))

# Build the design matrix
gse140719_design <- model.matrix(~gse140719_group)

# Fit the linear model and find DE (differentially expressed) probes
gse140719_fit <- lmFit(expr_mat_mimat, gse140719_design)
gse140719_fit <- eBayes(gse140719_fit, trend=TRUE)
summary(decideTests(gse140719_fit))
gse140719_top <- topTable(gse140719_fit, coef="gse140719_grouptumor", number=Inf, adjust.method="BH")

# Convert the result matrix into tibble
gse140719_result <- as_tibble(gse140719_top, rownames = "MIMAT_ID") %>%
    dplyr::select(MIMAT_ID, logFC, adj.P.Val)

# Get the significant IDs (adjusted p-value < 0.05 and fold change > 2)
gse140719_significant = gse140719_result %>%
    filter(abs(logFC) > 1 & adj.P.Val < 0.05)

#######################################
# Get overlapping significant MIMAT IDs
gse140719_significant_vec = pull(gse140719_significant, MIMAT_ID)
gse43796_significant_vec = pull(gse43796_significant, MIMAT_ID)
overlapping_MIMAT_ID = intersect(gse140719_significant_vec, gse43796_significant_vec)
print(overlapping_MIMAT_ID)
print(length(overlapping_MIMAT_ID))

# Check if overlapping IDs have same trend in both datasets
gse43796_11IDs = gse43796_significant %>%
    filter(MIMAT_ID %in% overlapping_MIMAT_ID)  %>%
    dplyr::rename(gse43796_logFC = logFC) %>%
    dplyr::select(MIMAT_ID, gse43796_logFC)

gse140719_11IDs = gse140719_significant %>%
    filter(MIMAT_ID %in% overlapping_MIMAT_ID)  %>%
    dplyr::rename(gse140719_logFC = logFC) %>%
    dplyr::select(MIMAT_ID, gse140719_logFC)

combined_11IDs = inner_join(gse43796_11IDs, gse140719_11IDs, join_by(MIMAT_ID))

#######################################
# Make a scatter plot for all the fold change for both data sets
gse43796_p_values = gse43796_result %>%
    mutate(gse43796_log_p = logFC) %>%
    select(MIMAT_ID, gse43796_log_p)

gse140719_p_values = gse140719_result %>%
    mutate(gse140719_log_p = logFC) %>%
    select(MIMAT_ID, gse140719_log_p)

merged_p_values <- inner_join(gse43796_p_values, gse140719_p_values, by = "MIMAT_ID")

cor.test(pull(merged_p_values, gse43796_log_p), pull(merged_p_values, gse140719_log_p))

ggplot(merged_p_values, aes(x = gse43796_log_p, y = gse140719_log_p)) +
    geom_point() +
    theme_bw()

#######################################
# Use multiMiR and miRBaseConverter to map miRNAs to their regulating mRNAs
# Install required library
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(c("multiMiR", "miRBaseConverter"))

library(multiMiR)
library(miRBaseConverter)

# Convert MIMAT IDs to mature names (e.g., hsa-miR-21-5p)
miRNA_names <- miRNA_AccessionToName(overlapping_MIMAT_ID, targetVersion = "v22")
miRNA_targets <- get_multimir(org = "hsa", mirna = miRNA_names$TargetName, table = "validated")
target_tibble <- as_tibble(miRNA_targets@data) %>%
    filter(support_type == "Functional MTI")

target_vec = pull(target_tibble, target_entrez)
print(length(unique(target_vec))) # Got 304 genes

########################################
# Compare targets with overlapping mRNA IDs and each data set IDs separately
# Load significant ID RDS object
emexp_1914_significant_vec <- readRDS("emexp_1914_significant_vec.rds") # 2540 genes
gse43795_significant_vec <- readRDS("gse43795_significant_vec.rds") # 4854 genes
overlapping_Entrez_ID <- readRDS("overlapping_Entrez_ID.rds") # 681 genes

target_vec = unique(target_vec)

# Get overlapping genes between emexp-1914 and targets, named overlapping_1
overlapping_1 = intersect(target_vec, emexp_1914_significant_vec)
print(overlapping_1)
print(length(overlapping_1)) # Got 53

# Get overlapping genes between gse43795 and targets, named overlapping_2
overlapping_2 = intersect(target_vec, gse43795_significant_vec)
print(overlapping_2)
print(length(overlapping_2)) # Got 114

# Get overlapping genes between overlapping_entrez_ID and targets, named overlapping_3
overlapping_3 = intersect(target_vec, overlapping_Entrez_ID)
print(overlapping_3)
print(length(overlapping_3)) # Got 18