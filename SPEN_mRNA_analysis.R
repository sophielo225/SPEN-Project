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

probeEnsemblRef <- getRefInfo(annotationPackagePrefix, "ENSEMBL")

# Map Ensembl gene ID to expression matrix
ensembl_ids <- probeEnsemblRef$ensembl_id[
    match(rownames(exprs_mat), probeEnsemblRef$probe_id)
]
# sum(is.na(ensembl_ids)) -> 3819
# This means that there are 3819 probes not annotated
# can we just ignore and remove them?

# Map ensembl IDs to expression matrix, and keep illumina probes when there's
# no ensembl genes annotated
rownames(exprs_mat) <- ifelse(
    is.na(ensembl_ids),
    rownames(exprs_mat),
    ensembl_ids
)

# Average duplicate probes through limma avereps()
gene_ids <- sub("\\.\\d+$", "", rownames(exprs_mat))
exprs_mat <- avereps(exprs_mat, ID = gene_ids)

# Select only rows that have ensembl gene IDs
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
    rownames_to_column(var = "Ensembl_ID") %>%
    as_tibble() %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    dplyr::select(Ensembl_ID, logFC, adj.P.Val)

########################################
# E-MEXP-1914
emexp_1914 = readRDS("E-MEXP-1914.eSet.rds")

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

# Build the design matrix
design <- model.matrix(~ 1, data = data.frame(array = colnames(MA)))

# Fit the model and do differential expression analysis
fit <- lmFit(MA, design)
fit <- eBayes(fit)
emexp_1914_results <- topTable(fit, coef = 1, number = Inf, adjust.method="BH")

# Filter the significant genes
emexp_1914_significant <- emexp_1914_results %>%
    as_tibble() %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    filter(grepl("^ENST", `Composite.Element.Database.Entry.ensembl.`))
########################################
# old script (manually constructed)
# Download idf, sdrf, and raw files
dir.create("E-MEXP-1914", showWarnings = FALSE)
dest_dir <- "E-MEXP-1914"

file_urls <- c(
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/E-MEXP-1914.idf.txt",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/E-MEXP-1914.sdrf.txt",
    "https://ftp.ebi.ac.uk/biostudies/fire/A-MEXP-/703/A-MEXP-703/Files/A-MEXP-703.adf.txt",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/Agilent46820-TSPP2.gpr",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/Agilent46821-TSPP3.gpr",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/Agilent46823-TSPP6.gpr",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/Agilent46754-TSPP1.gpr",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/914/E-MEXP-1914/Files/Agilent46822-TSPP5.gpr"
)

for (url in file_urls) {
    dest_file <- file.path(dest_dir, basename(url))
    download.file(url, destfile = dest_file, mode = "wb")
    message("Downloaded: ", dest_file)
}

# Read metadata and feature data
e_mexp_1914_metadata = read_tsv("E-MEXP-1914/E-MEXP-1914.sdrf.txt")
e_mexp_1914_idf = read_tsv("E-MEXP-1914/E-MEXP-1914.idf.txt")

adf_file <- "E-MEXP-1914/A-MEXP-703.adf.txt"

# Find where the [main] section begins
lines <- readLines(adf_file)
start <- grep("^\\[main\\]", lines)

# Read the feature table starting *after* the [main] line
e_mexp_1914_feature_data <- read.delim(adf_file, skip = start, header = TRUE, sep="\t", quote="")

# Generate target file that contains the info for reading the files
# --- Define the function to read a single GPR file ---
read_custom_gpr_v2 <- function(filename) {
    
    # Use read.table with robust settings for tab-delimited files
    data <- read.table(
        filename, 
        header = TRUE, 
        sep = "\t", 
        skip = 0,               # Set to 0 since you believe the header is on line 1
        quote = "",             # Ignore quotes within fields
        comment.char = "",      # CRUCIAL: Treats all lines as data, ignoring possible comment characters
        stringsAsFactors = FALSE 
    )
    
    # Note: Column names will have periods instead of spaces/colons 
    # (e.g., GenePix:F635 Mean becomes GenePix.F635.Mean)
    
    # 2. Extract the required intensity data columns
    R <- data$GenePix.F635.Mean    # Red/Cy5 Foreground Mean
    G <- data$GenePix.F532.Mean    # Green/Cy3 Foreground Mean
    Rb <- data$GenePix.B635.Mean   # Red/Cy5 Background Mean
    Gb <- data$GenePix.B532.Mean   # Green/Cy3 Background Mean
    
    # 3. Extract metadata columns
    ID <- data$CompositeSequence.Identifier
    Flags <- data$GenePix.Flags
    
    return(list(
        R = R, G = G, Rb = Rb, Gb = Gb, ID = ID, Flags = Flags
    ))
}

# --- Set the variables and run the analysis ---
gpr_folder <- "E-MEXP-1914/" 

# 1. Get the list of file paths
targets <- list.files(
    path = gpr_folder, 
    pattern = "\\.gpr$", 
    full.names = TRUE
)

# 2. Apply the custom function to all files
data_list <- lapply(targets, read_custom_gpr_v2)

# 3. Combine the list elements into a single RGList object
RG <- new("RGList")
RG$R <- sapply(data_list, function(x) x$R)
RG$G <- sapply(data_list, function(x) x$G)
RG$Rb <- sapply(data_list, function(x) x$Rb)
RG$Gb <- sapply(data_list, function(x) x$Gb)

# Add probe metadata (using the ID and Flags from the first file)
RG$genes <- data.frame(
    ID = data_list[[1]]$ID, 
    Flags = data_list[[1]]$Flags 
)

# Add file names
RG$targets <- data.frame(FileName = targets)

# Check the structure
print(RG)

# Background Correction (using "normexp" is recommended)
MA <- backgroundCorrect(RG, method = "normexp")

# Define a quality weight vector (1 = keep, 0 = remove)
weights <- as.numeric(RG$genes$Flags >= 0)

# Optional: Add the weights to the MA object
MA$weights <- weights

# Normalize the data
MA.norm <- normalizeWithinArrays(MA, method = "loess", weights = MA$weights)

# Fit the linear model to the normalized data (MA.norm)
fit <- lmFit(MA.norm)

# Apply Empirical Bayes Moderation
fit <- eBayes(fit)

# Extract the results
results <- topTable(
    fit, 
    number = Inf, 
    adjust.method = "fdr",
    # The coefficient of interest is the intercept, which represents the 
    # log2(SPEN/Benign) comparison. By default, limma uses column 1.
    coef = 1 
)

# View the top differentially expressed genes
head(results)