# Load all the required libraries
library(GEOquery)
library(tidyverse)
library(limma)
library(oligo)
library(R.utils)
library(pd.mirna.4.0)

# Process GSE43796 data set
# Get metadata and feature data
gse43796 <- getGEO("GSE43796", GSEMatrix = TRUE)

gse43796_metadata <- pData(gse43796[[1]]) %>%
    as_tibble(rownames = "Sample_ID") %>%
    unite(Sample_ID, Sample_ID, description, sep="_") %>%
    filter(source_name_ch1 %in% c("solid-pseudopapillary neoplasm", "non-neoplastic pancreas")) 

gse43796_feature_data <- fData(gse43796[[1]]) %>%
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

# Perform background correction
RG_bc <- limma::backgroundCorrect(RG, method="normexp")  

# Perform normalization between arrays
MA <- normalizeBetweenArrays(RG_bc, method="quantile")

# Remove directory names from the column names in the expression matrix
colnames(MA$E) <- basename(colnames(MA$E))

# Only keep the samples we care about
MA <- MA[,pull(gse43796_metadata, Sample_ID)]

# Only keep regular probes (remove control probes)
MA <- MA[MA$genes$ControlType==0, ]

# Map MIMAT IDs to each probe name in the expression matrix
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

# Fit the linear model and find differentially expressed probes
fit <- lmFit(MA_collapse, design)
fit <- eBayes(fit, trend=TRUE)
summary(decideTests(fit))
top <- topTable(fit, coef="grouptumor", number=Inf, adjust.method="BH")

# Convert the result matrix into tibble
gse43796_result <- as_tibble(top) %>%
    dplyr::select(MIMAT_ID, SystematicName, logFC, adj.P.Val)

# Get the significant IDs (adjusted p-value < 0.05 and log fold change > 1)
gse43796_significant = gse43796_result %>%
    filter(abs(logFC) > 1 & adj.P.Val < 0.05)

#######################################
# Process GSE140719 data set
# Get metadata and feature data
gse140719 <- getGEO("GSE140719", GSEMatrix = TRUE)

gse140719_metadata <- pData(gse140719[[1]]) %>%
    as_tibble(rownames = "Sample_ID") %>%
    filter(grepl("localized|Normal", title)) %>%
    dplyr::select(Sample_ID, title)

gse140719_feature_data <- fData(gse140719[[1]]) %>%
    as_tibble(rownames = "ID_REF") %>%
    dplyr::select(ID_REF, Accession)

# Get CEL files
gse_id <- "GSE140719"
cel_dir <- file.path(gse_id, "CEL")

# Unzip all .cel.gz files
cel_archives <- list.files(cel_dir, pattern = "\\.cel\\.gz$", full.names = TRUE)
for (f in cel_archives) {
    R.utils::gunzip(f, overwrite = TRUE, remove = TRUE)
}

cel_files <- list.celfiles(cel_dir, full.names = TRUE)
stopifnot(length(cel_files) > 0)

# Read raw CELs without annotation
raw <- read.celfiles(cel_files, pkgname = "pd.mirna.4.0")

# Background correction + quantile normalization + summarization
eset <- rma(raw)                  
expr_mat <- Biobase::exprs(eset)

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

# Fit the linear model and find differentially expressed probes
gse140719_fit <- lmFit(expr_mat_mimat, gse140719_design)
gse140719_fit <- eBayes(gse140719_fit, trend=TRUE)
summary(decideTests(gse140719_fit))
gse140719_top <- topTable(gse140719_fit, coef="gse140719_grouptumor", number=Inf, adjust.method="BH")

# Convert the result matrix into tibble
gse140719_result <- as_tibble(gse140719_top, rownames = "MIMAT_ID") %>%
    dplyr::select(MIMAT_ID, logFC, adj.P.Val)

# Get the significant IDs (adjusted p-value < 0.05 and log fold change > 1)
gse140719_significant <- gse140719_result %>%
    filter(abs(logFC) > 1 & adj.P.Val < 0.05)

#######################################
# Get overlapping significant MIMAT IDs
gse140719_significant_vec <- pull(gse140719_significant, MIMAT_ID)
gse43796_significant_vec <- pull(gse43796_significant, MIMAT_ID)
overlapping_MIMAT_ID <- intersect(gse140719_significant_vec, gse43796_significant_vec)
print(overlapping_MIMAT_ID)
print(length(overlapping_MIMAT_ID)) # 11 overlapping miRNAs

# Check if overlapping IDs have same trend of fold change in both data sets
gse43796_11IDs <- gse43796_significant %>%
    filter(MIMAT_ID %in% overlapping_MIMAT_ID)  %>%
    dplyr::rename(gse43796_logFC = logFC) %>%
    dplyr::select(MIMAT_ID, gse43796_logFC)

gse140719_11IDs <- gse140719_significant %>%
    filter(MIMAT_ID %in% overlapping_MIMAT_ID)  %>%
    dplyr::rename(gse140719_logFC = logFC) %>%
    dplyr::select(MIMAT_ID, gse140719_logFC)

combined_11IDs <- inner_join(gse43796_11IDs, gse140719_11IDs, join_by(MIMAT_ID))

# Divide 11 miRNAs into two groups based on the sign of the fold change 
positive_FC_IDs <- combined_11IDs %>%
    filter(gse43796_logFC > 0 & gse140719_logFC > 0) %>%
    pull(MIMAT_ID) # Up-regulated miRNAs

negative_FC_IDs <- combined_11IDs %>%
    filter(gse43796_logFC < 0 & gse140719_logFC < 0) %>%
    pull(MIMAT_ID) # Down-regulated miRNAs

# Save RDS object to be used in other files
saveRDS(overlapping_MIMAT_ID, "overlapping_MIMAT_ID.rds")
saveRDS(positive_FC_IDs, "positive_FC_ID.rds")
saveRDS(negative_FC_IDs, "negative_FC_ID.rds")