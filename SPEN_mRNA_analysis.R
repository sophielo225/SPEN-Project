# This file analyzes two mRNA sequencing data sets(GSE43795, E-MEXP-1914)
# GSE43795 is from GEO, and E-MEXP-1914 is from ArrayExpress

# Install needed packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")

# Load libraries
library(GEOquery)
library(tidyverse)
library(limma)
library(R.utils)

# GSE43795
# Get metadata and feature data
gse43795 = getGEO("GSE43795", GSEMatrix = TRUE)

gse43795_metadata = pData(gse43795[[1]]) %>%
    as_tibble(rownames = "Sample_ID") %>%
    filter(source_name_ch1 == "solid-pseudopapillary neoplasm" | source_name_ch1 == "non-neoplastic pancreas")

gse43795_feature_data = fData(gse43795[[1]]) %>%
    as_tibble(rownames = "ID_REF") %>%
    dplyr::select(ID, Entrez_Gene_ID, Symbol)

# Download all supplementary files for the series (includes the .tar raw file)
getGEOSuppFiles("GSE43795")

# Set the file path
tarfile <- "GSE43795/GSE43795_RAW.tar"

# Unpack the tar file into a subfolder
untar(tarfile, exdir = "GSE43795/raw_data")

# Gunzip them if they are compressed
txt_files <- list.files("GSE43795/raw_data", pattern = "\\.gz$", full.names = TRUE)
sapply(txt_files, R.utils::gunzip, overwrite=TRUE)

raw_files <- list.files("GSE43795/raw_data", pattern = "\\.txt$", full.names = TRUE)

# Path to the file
gzfile <- "GSE43795/GSE43795_non_normalized.txt.gz"

# Unzip into the same folder
gunzip(gzfile, overwrite = TRUE)

# Read the non-normalized TXT file
data <- read.delim("GSE43795/GSE43795_non_normalized.txt", stringsAsFactors = FALSE)

# Expression columns: S.*, A.*, E.*, N.*
expr_cols <- grep("^[SANE]\\.", colnames(data))

# Detection P-values: Detection.Pval.*
det_cols <- grep("^Detection.Pval", colnames(data))

# Create expression and detection matrices
exprs_mat <- as.matrix(data[, expr_cols])
detP_mat <- as.matrix(data[, det_cols])
rownames(exprs_mat) <- data$ID_REF
rownames(detP_mat) <- data$ID_REF
colnames(detP_mat) <- colnames(exprs_mat)

# Perform log2 transformation, background correction, quantile normalization 
exprs_mat <- neqc(x = exprs_mat, detection.p = detP_mat)

# Only keep the samples we care about
cols_to_keep <- c(paste0("S.", 1:13), paste0("N.", 1:5))
exprs_mat <- exprs_mat[, cols_to_keep]
detP_mat <- detP_mat[, cols_to_keep]

# Exclude expression rows where all detection p-values are < 0.05
expressed <- apply(detP_mat < 0.05, 1, any)
exprs_mat <- exprs_mat[expressed,]

# BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)

# TODO: Move to a better spot?
getRefInfo <- function(annotationPackagePrefix, suffix, secondColumnName=NULL) {
    info <- get(paste0(annotationPackagePrefix, suffix))
    info <- info[mappedkeys(info)] # Returns the subset of mapped keys.
    
    df <- as.data.frame(info)
    
    if (!is.null(secondColumnName)) {
        colnames(df) <- c("Illumina_ID", secondColumnName)
    }
    
    return(df)
}

annotationPackagePrefix <- "illuminaHumanv4"

probeQualityRef <- getRefInfo(annotationPackagePrefix, "PROBEQUALITY")
probeReporterRef <- getRefInfo(annotationPackagePrefix, "REPORTERGROUPNAME", "Probe_Reporter_Type")

controlProbes <- probeReporterRef$Illumina_ID[which(probeReporterRef$Probe_Reporter_Type == "negative")]
controlProbes <- intersect(controlProbes, rownames(exprs_mat))

perfectProbes <- probeQualityRef$IlluminaID[which(grepl("Perfect", probeQualityRef$ProbeQuality))]
goodProbes <- probeQualityRef$IlluminaID[which(grepl("Good", probeQualityRef$ProbeQuality))]
effectiveProbes <- c(perfectProbes, goodProbes)

signalExprProbes <- intersect(rownames(exprs_mat), effectiveProbes)
signalExprProbes <- setdiff(signalExprProbes, controlProbes)

exprs_mat <- exprs_mat[signalExprProbes,]
detP_mat <- detP_mat[signalExprProbes,]

# Depends on the gene IDs we want to use, I've included both Entrez and Ensembl
probeEnsemblRef <- getRefInfo(annotationPackagePrefix, "ENSEMBL")
probeEntrezRef <- getRefInfo(annotationPackagePrefix, "ENTREZID")

# Map Ensembl gene ID to expression matrix
ensembl_ids <- probeEnsemblRef$ensembl_id[
    match(rownames(exprs_mat), probeEnsemblRef$probe_id)
]

# Map ensembl IDs to expression matrix, and keep illumina probes when there's
# no ensembl genes annotated
rownames(exprs_mat) <- ifelse(
    is.na(ensembl_ids),
    rownames(exprs_mat),
    ensembl_ids
)

# Map Entrez IDs to expression matrix
entrez_ids <- probeEntrezRef$gene_id[
    match(rownames(exprs_mat), probeEntrezRef$probe_id)
]

# Map Entrez IDs to expression matrix, and keep illumina probes when there's
# no Entrez genes annotated
rownames(exprs_mat) <- ifelse(
    is.na(entrez_ids),
    rownames(exprs_mat),
    entrez_ids
)

# Average duplicate probes through limma avereps()
gene_ids <- sub("\\.\\d+$", "", rownames(exprs_mat))
exprs_mat <- avereps(exprs_mat, ID = gene_ids)

# Select only rows that have ensembl gene IDs
exprs_mat <- exprs_mat[!grepl("^ILMN_", rownames(exprs_mat)), ]

# Select only rows that have Entrez gene IDs
exprs_mat <- exprs_mat[!grepl("^ILMN_", rownames(exprs_mat)), ]

# Create the sample group factors
group <- ifelse(grepl("^S", colnames(exprs_mat)), "tumor", "normal")
group <- factor(group)

# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit the linear model
fit <- lmFit(exprs_mat, design)

# Perform differential expression analysis and get the number of differentially expressed probes
contrast_matrix <- makeContrasts(tumor - normal, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))
gse43795_results <- topTable(fit2, number=Inf, adjust.method="BH")

# Keep only significant rows (adj.p < 0.05 and fold change > 1)
gse43795_significant = gse43795_results %>%
    rownames_to_column(var = "Entrez_ID") %>%
    as_tibble() %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    dplyr::select(Entrez_ID, logFC, adj.P.Val)

########################################
# TODO: Maybe move this to somewhere else
BiocManager::install("hgug4110b.db")
library(hgug4110b.db)

# E-MEXP-1914
emexp_1914 = readRDS("E-MEXP-1914.eSet.rds")

## Bimap interface:
x <- hgug4110bENTREZID
# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(x)
# Convert it to a list
xx <- as.list(x[mapped_probes])
# Convert it to a tibble
emexp_1914_probe_annotation <- enframe(xx, name = "probe_ID", value = "Entrez_ID") %>%
    unnest(Entrez_ID)

# Manually construct it into RGList
RG <- new("RGList", list(
    R  = assayData(emexp_1914)$R,
    G  = assayData(emexp_1914)$G,
    Rb = assayData(emexp_1914)$Rb,
    Gb = assayData(emexp_1914)$Gb
))

RG$genes <- fData(emexp_1914)
RG$targets <- pData(emexp_1914)

# Background correction
RG <- backgroundCorrect(RG, method = "normexp", offset = 50)

# Within-array normalization
MA <- normalizeWithinArrays(RG, method = "loess")

# Between-array normalization
MA <- normalizeBetweenArrays(MA, method = "Aquantile")

# Map Entrez gene IDs to feature data
feature_tibble <- rownames_to_column(emexp_1914@featureData@data, var = "ref_num") %>%
    select(ref_num, `Reporter.Database.Entry.agilent_probe.`)

# Map Entrez gene IDs to expression matrix
expression_tibble <- as.data.frame(MA$M) %>%
    rownames_to_column(var = "ref_num")

# Join expression tibble with probe annotation tibble
combined_tibble <- inner_join(expression_tibble, feature_tibble, join_by(ref_num)) %>%
    select(`Reporter.Database.Entry.agilent_probe.`, `Agilent46822-TSPP5`:`Agilent46820-TSPP2`) %>%
    rename(probe_ID = `Reporter.Database.Entry.agilent_probe.`) %>%
    filter(grepl("^A_", probe_ID))

# Map probe-mapped expression tibble with Entrez gene ID
annotated_tibble <- inner_join(combined_tibble, emexp_1914_probe_annotation, join_by(probe_ID)) %>%
    select(Entrez_ID, `Agilent46822-TSPP5`:`Agilent46820-TSPP2`)

# Since matrix does not allow duplicate row names, so I exclude Entrez_ID and 
# convert the tibble back to matrix. In this case, I can do avereps(), since 
# avereps() does not require row names for the expression matrix
annotated_expr_mat <- as.matrix(annotated_tibble %>% select(-Entrez_ID))
expr_gene <- avereps(annotated_expr_mat, ID = annotated_tibble$Entrez_ID)

# Build the design matrix
design <- model.matrix(~ 1, data = data.frame(array = colnames(expr_gene)))

# Fit the model and do differential expression analysis
fit <- lmFit(expr_gene, design)
fit <- eBayes(fit)
emexp_1914_results <- topTable(fit, coef = 1, number = Inf, adjust.method="BH")

# Filter the significant genes
emexp_1914_significant <- emexp_1914_results %>%
    rownames_to_column(var = "Entrez_ID") %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Get overlapping gene IDs from both datasets
emexp_1914_significant_vec = pull(emexp_1914_significant, Entrez_ID)
gse43795_significant_vec = pull(gse43795_significant, Entrez_ID)
overlapping_Entrez_ID = intersect(emexp_1914_significant_vec, gse43795_significant_vec)
print(overlapping_Entrez_ID)
print(length(overlapping_Entrez_ID))

###