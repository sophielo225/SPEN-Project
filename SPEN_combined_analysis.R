# Install required library
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(c("multiMiR", "miRBaseConverter"))

# Use multiMiR and miRBaseConverter to map miRNAs to their regulating mRNAs
library(multiMiR)
library(miRBaseConverter)

# Load significant ID RDS object
emexp_1914_significant_vec <- readRDS("emexp_1914_significant_vec.rds") # 2540 genes
gse43795_significant_vec <- readRDS("gse43795_significant_vec.rds") # 4854 genes
overlapping_Entrez_ID <- readRDS("overlapping_Entrez_ID.rds") # 681 genes
same_trend_overlapping_Entrez_ID <- readRDS("same_trend_overlapping_Entrez_ID.rds") # 336 genes
overlapping_MIMAT_ID <- readRDS("overlapping_MIMAT_ID.rds") # 11 miRNAs

# Convert MIMAT IDs to mature names (e.g., hsa-miR-21-5p)
miRNA_names <- miRNA_AccessionToName(overlapping_MIMAT_ID, targetVersion = "v22")
miRNA_targets <- get_multimir(org = "hsa", mirna = miRNA_names$TargetName, table = "validated")
target_tibble <- as_tibble(miRNA_targets@data) %>%
    filter(support_type == "Functional MTI")

target_vec <- pull(target_tibble, target_entrez)
print(length(unique(target_vec))) # Got 304 genes

# Compare targets with overlapping mRNA IDs and each data set IDs separately
target_vec <- unique(target_vec)

# Get overlapping genes between emexp-1914 and targets, named overlapping_1
overlapping_1 <- intersect(target_vec, emexp_1914_significant_vec)
print(overlapping_1)
print(length(overlapping_1)) # Got 53

# Get overlapping genes between gse43795 and targets, named overlapping_2
overlapping_2 <- intersect(target_vec, gse43795_significant_vec)
print(overlapping_2)
print(length(overlapping_2)) # Got 114

# Get overlapping genes between overlapping_entrez_ID and targets, named overlapping_3
overlapping_3 <- intersect(target_vec, overlapping_Entrez_ID)
print(overlapping_3)
print(length(overlapping_3)) # Got 18

# Get overlapping genes between same_trend_overlapping_Entrez_ID and targets, named overlapping_4
overlapping_4 <- intersect(target_vec, same_trend_overlapping_Entrez_ID)
print(overlapping_4)
print(length(overlapping_4)) # Got 10