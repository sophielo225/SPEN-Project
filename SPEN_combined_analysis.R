# Use multiMiR and miRBaseConverter to map miRNAs to their regulating mRNAs
library(multiMiR)
library(miRBaseConverter)
library(tidyverse)

# Load significant ID RDS object
emexp_1914_significant_vec <- readRDS("emexp_1914_significant_vec.rds") # 2540 genes
gse43795_significant_vec <- readRDS("gse43795_significant_vec.rds") # 4854 genes
same_trend_overlapping_Entrez_ID <- readRDS("same_trend_overlapping_Entrez_ID.rds") # 336 genes
overlapping_MIMAT_ID <- readRDS("overlapping_MIMAT_ID.rds") # 11 miRNAs

# Convert MIMAT IDs to mature names (e.g., hsa-miR-21-5p) using miRBaseConverter
miRNA_names <- miRNA_AccessionToName(overlapping_MIMAT_ID, targetVersion = "v22")

# Use multiMiR to map miRNAs to their validated target mRNAs
miRNA_targets <- get_multimir(org = "hsa", mirna = miRNA_names$TargetName, table = "validated")
target_tibble <- as_tibble(miRNA_targets@data) %>%
    filter(support_type == "Functional MTI")

target_genes = dplyr::select(target_tibble, target_symbol, target_entrez) %>%
  distinct()
target_vec <- pull(target_genes, target_entrez)
print(length(unique(target_vec))) # Got 304 target genes

# Compare targets with overlapping mRNA IDs and each data set IDs separately
target_vec <- unique(target_vec) # n = 304

# Get overlapping genes between emexp-1914 and targets
overlapping_1 <- intersect(target_vec, emexp_1914_significant_vec)
print(overlapping_1)
print(length(overlapping_1)) # Got 53 (out of 2,540)

# Get overlapping genes between gse43795 and targets
overlapping_2 <- intersect(target_vec, gse43795_significant_vec)
print(overlapping_2)
print(length(overlapping_2)) # Got 114 (out of 4,854)
# These 114 genes' Entrez IDs were passed in to do gene-set enrichment analysis
# The resulting tsv file is called "GSEA_result_114_genes.tsv" and saved under "pathway_results" folder

# Get overlapping genes between same_trend_overlapping_Entrez_ID and targets
overlapping_3 <- intersect(target_vec, same_trend_overlapping_Entrez_ID)
print(overlapping_3)
print(length(overlapping_3)) # Got 10 (out of 681)
# These 10 genes' Entrez IDs were passed in to do gene-set enrichment analysis
# The resulting tsv file is called "GSEA_result_10_genes.tsv" and saved under "pathway_results" folder

########################################
# Expected directionality analysis
# Compare up- or down-regulated miRNA IDs with up- or down- regulation of their target mRNA
# If a miRNA is up-regulated, then its target mRNA should be down-regulated
# If a miRNA is down-regulated, then its target mRNA should be up-regulated

# Load RDS object
positive_FC_IDs <- readRDS("positive_FC_ID.rds")
negative_FC_IDs <- readRDS("negative_FC_ID.rds")
gse43795_significant_vec <- readRDS("gse43795_significant_vec.rds")
gse43795_significant <- readRDS("gse43795_significant_tibble.rds")
emexp_1914_significant <- readRDS("emexp_1914_significant_tibble.rds")
overlapping_336_genes <- readRDS("overlapping_336_genes.rds")

# Convert Convert MIMAT IDs to mature names
positive_FC_miRNA <- miRNA_names %>%
    filter(Accession %in% positive_FC_IDs) %>%
    pull(TargetName)

negative_FC_miRNA <- miRNA_names %>%
    filter(Accession %in% negative_FC_IDs) %>%
    pull(TargetName)

gse43795_positive_FC_gene <- gse43795_significant %>%
  filter(logFC > 1) %>%
  pull(Entrez_Gene_ID) # 2754

gse43795_negative_FC_gene <- gse43795_significant %>%
  filter(logFC < -1) %>%
  pull(Entrez_Gene_ID) # 2100

emexp_1914_positive_FC_gene <- emexp_1914_significant %>%
  filter(logFC > 1) %>%
  pull(Entrez_Gene_ID) # 1419

emexp_1914_negative_FC_gene <- emexp_1914_significant %>%
  filter(logFC < -1) %>%
  pull(Entrez_Gene_ID) # 1121

both_mRNA_positive_FC_gene <- intersect(gse43795_positive_FC_gene, emexp_1914_positive_FC_gene) # 187

both_mRNA_negative_FC_gene <- intersect(gse43795_negative_FC_gene, emexp_1914_negative_FC_gene) # 149

miRNA_targets <- get_multimir(org = "hsa", mirna = c(positive_FC_miRNA, negative_FC_miRNA), table = "validated")
miRNA_targets <- as_tibble(miRNA_targets@data) %>%
  filter(support_type == "Functional MTI") %>%
  distinct(mature_mirna_id, target_symbol, target_entrez) # 351

positive_miRNA_targets <- filter(miRNA_targets, mature_mirna_id %in% positive_FC_miRNA) %>%
  pull(target_entrez) %>%
  unique() # 191 target genes that should be, in theory, but they may not be in our differentally expressed gene list.

negative_miRNA_targets <- filter(miRNA_targets, mature_mirna_id %in% negative_FC_miRNA) %>%
  pull(target_entrez) %>%
  unique() # 132 target genes that should, in theory, be downregulated, but they may not be in our differentially expressed gene list.

intersect(gse43795_positive_FC_gene, c(positive_miRNA_targets, negative_miRNA_targets)) %>%
  length() # 59 all relationships
intersect(gse43795_positive_FC_gene, negative_miRNA_targets) %>%
  length() # 34 positive inverse/expected
intersect(gse43795_positive_FC_gene, positive_miRNA_targets) %>%
  length() # 28 positive same
intersect(gse43795_negative_FC_gene, positive_miRNA_targets) %>%
  length() # 36 negative inverse/expected
intersect(gse43795_negative_FC_gene, negative_miRNA_targets) %>%
  length() # 22 negative same

intersect(both_mRNA_positive_FC_gene, c(positive_miRNA_targets, negative_miRNA_targets)) %>%
  length() # 7 all relationships
intersect(both_mRNA_positive_FC_gene, negative_miRNA_targets) %>%
  length() # 5 positive inverse/expected
intersect(both_mRNA_positive_FC_gene, positive_miRNA_targets) %>%
  length() # 2 positive same
intersect(both_mRNA_negative_FC_gene, positive_miRNA_targets) %>%
  length() # 1 negative inverse/expected
intersect(both_mRNA_negative_FC_gene, negative_miRNA_targets) %>%
  length() # 2 negative same

gse43795_positive_inverse_genes <- intersect(gse43795_positive_FC_gene, negative_miRNA_targets)
print(gse43795_positive_inverse_genes)

# These 36 genes' Entrez IDs were passed in to do gene-set enrichment analysis
# The resulting tsv file is called "gse43795_negative_FC_GSEA_result.tsv" and saved under "pathway_results" folder

gse43795_negative_inverse_genes <- intersect(gse43795_negative_FC_gene, positive_miRNA_targets)
print(gse43795_negative_inverse_genes)

# These 34 genes' Entrez IDs were passed in to do gene-set enrichment analysis
# The resulting tsv file is called "gse43795_positive_FC_GSEA_result.tsv" and saved under "pathway_results" folder

both_mRNA_positive_inverse_genes <- intersect(both_mRNA_positive_FC_gene, negative_miRNA_targets)
print(both_mRNA_positive_inverse_genes)

# These 5 genes' Entrez IDs were passed in to do gene-set enrichment analysis
# But there's no identified pathways overlapping across these genes

both_mRNA_negative_inverse_genes <- intersect(both_mRNA_negative_FC_gene, positive_miRNA_targets)
print(both_mRNA_negative_inverse_genes)

# There's only one gene, so we didn't do a gene-set enrichment analysis