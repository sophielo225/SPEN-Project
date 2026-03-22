# Use multiMiR and miRBaseConverter to map miRNAs to their regulating mRNAs
library(multiMiR)
library(miRBaseConverter)
library(tidyverse)

# Load significant ID RDS object
emexp_1914_significant_vec <- readRDS("emexp_1914_significant_vec.rds") # 2540 genes
gse43795_significant_vec <- readRDS("gse43795_significant_vec.rds") # 4854 genes
overlapping_Entrez_ID <- readRDS("overlapping_Entrez_ID.rds") # 681 genes
same_trend_overlapping_Entrez_ID <- readRDS("same_trend_overlapping_Entrez_ID.rds") # 336 genes
overlapping_MIMAT_ID <- readRDS("overlapping_MIMAT_ID.rds") # 11 miRNAs

# Convert MIMAT IDs to mature names (e.g., hsa-miR-21-5p) using miRBaseConverter
miRNA_names <- miRNA_AccessionToName(overlapping_MIMAT_ID, targetVersion = "v22")

# Use multiMiR to map miRNAs to their validated target mRNAs
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

########################################
# Compare up- or down-regulated miRNA IDs with up- or down- regulation of their target mRNA
# If a miRNA is up-regulated, then its target mRNA should be down-regulated
# If a miRNA is down-regulated, then its target mRNA should be up-regulated

# Load RDS object
positive_FC_IDs <- readRDS("positive_FC_ID.rds")
negative_FC_IDs <- readRDS("negative_FC_ID.rds")
gse43795_significant_vec <- readRDS("gse43795_significant_vec.rds")
gse43795_significant <- readRDS("gse43795_significant_tibble.rds")

# Convert Convert MIMAT IDs to mature names
positive_FC_miRNA <- miRNA_names %>%
    filter(Accession %in% positive_FC_IDs) %>%
    pull(TargetName)

negative_FC_miRNA <- miRNA_names %>%
    filter(Accession %in% negative_FC_IDs) %>%
    pull(TargetName)

# Get each group's target mRNA
positive_miRNA_targets <- get_multimir(org = "hsa", mirna = positive_FC_miRNA, table = "validated")
positive_target_tibble <- as_tibble(positive_miRNA_targets@data) %>%
    filter(support_type == "Functional MTI")

negative_miRNA_targets <- get_multimir(org = "hsa", mirna = negative_FC_miRNA, table = "validated")
negative_target_tibble <- as_tibble(negative_miRNA_targets@data) %>%
    filter(support_type == "Functional MTI")

negative_target_vec <- unique(pull(positive_target_tibble, target_entrez))
print(length(negative_target_vec)) # Got 191 genes that should be down-regulated

positive_target_vec <- unique(pull(negative_target_tibble, target_entrez))
print(length(positive_target_vec)) # Got 132 genes that should be up-regulated

# Compare those target mRNAs to GSE43795 dataset to see the regulation of their fold change
negative_overlapping_gene <- intersect(negative_target_vec, gse43795_significant_vec)
positive_overlapping_gene <- intersect(positive_target_vec, gse43795_significant_vec)

negative_FC_gene <- gse43795_significant %>%
    filter(Entrez_ID %in% negative_overlapping_gene & logFC < 0) %>%
    pull(Entrez_ID)
print(length(negative_FC_gene)) # Got 36 actual target genes that are down-regulated

positive_FC_gene <- gse43795_significant %>%
    filter(Entrez_ID %in% positive_overlapping_gene & logFC > 0) %>%
    pull(Entrez_ID)
print(length(positive_FC_gene)) # Got 34 actual target genes that are up-regulated