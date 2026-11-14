# Use multiMiR and miRBaseConverter to map miRNAs to their regulating mRNAs
library(multiMiR)
library(miRBaseConverter)
library(tidyverse)

gene_annotations = mutate(gene_annotations, Entrez_Gene_ID = as.character(Entrez_Gene_ID))

# Convert MIMAT IDs to mature names (e.g., hsa-miR-21-5p) using miRBaseConverter
miRNA_names <- miRNA_AccessionToName(overlapping_MIMAT_ID, targetVersion = "v22")

# Use multiMiR to map miRNAs to their validated target mRNAs
miRNA_targets <- get_multimir(org = "hsa", mirna = miRNA_names$TargetName, table = "validated")
miRNA_targets <- as_tibble(miRNA_targets@data) %>%
  filter(support_type == "Functional MTI") %>%
  distinct(mature_mirna_id, target_entrez)

# Expected directionality analysis
# If a miRNA is up-regulated, then its target mRNA should be down-regulated
# If a miRNA is down-regulated, then its target mRNA should be up-regulated

# Convert Convert MIMAT IDs to mature names
positive_FC_miRNA <- miRNA_names %>%
    filter(Accession %in% positive_FC_IDs) %>%
    pull(TargetName) # 6

negative_FC_miRNA <- miRNA_names %>%
    filter(Accession %in% negative_FC_IDs) %>%
    pull(TargetName) # 5

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

positive_miRNA_targets_tibble <- filter(miRNA_targets, mature_mirna_id %in% positive_FC_miRNA) %>%
  filter(target_entrez != "")

positive_miRNA_targets <- pull(positive_miRNA_targets_tibble, target_entrez) %>%
  unique() # 190 target genes for 6 miRNA that were overexpressed in tumors.

negative_miRNA_targets_tibble <- filter(miRNA_targets, mature_mirna_id %in% negative_FC_miRNA) %>%
  filter(target_entrez != "")

negative_miRNA_targets <- pull(negative_miRNA_targets_tibble, target_entrez) %>%
  unique() # 131 target genes for 5 miRNA that were underexpressed in tumors.

intersect(c(positive_miRNA_targets, negative_miRNA_targets), c(both_mRNA_positive_FC_gene, both_mRNA_negative_FC_gene)) %>%
  length() # 10 all relationships
intersect(negative_miRNA_targets, both_mRNA_positive_FC_gene) %>%
  length() # 5 negative inverse/expected
intersect(negative_miRNA_targets, both_mRNA_negative_FC_gene) %>%
  length() # 2 negative same
intersect(positive_miRNA_targets, both_mRNA_negative_FC_gene) %>%
  length() # 1 positive inverse/expected
intersect(positive_miRNA_targets, both_mRNA_positive_FC_gene) %>%
  length() # 2 positive same

both_mRNA_negative_inverse_genes <- intersect(negative_miRNA_targets, both_mRNA_positive_FC_gene)

both_mRNA_positive_inverse_genes <- intersect(positive_miRNA_targets, both_mRNA_negative_FC_gene)

miRNA_targets %>%
  inner_join(gene_annotations, by = c("target_entrez" = "Entrez_Gene_ID")) %>%
  dplyr::rename(miRNA_ID = mature_mirna_id) %>%
  dplyr::rename(Target_Entrez_Gene_ID = target_entrez) %>%
  dplyr::rename(Target_Genome_Symbol = HGNC_Symbol) %>%
  mutate(miRNA_Expression_Direction = if_else(Target_Entrez_Gene_ID %in% positive_miRNA_targets, "positive", "negative")) %>%
  mutate(mRNA_Expression_Direction = if_else(Target_Entrez_Gene_ID %in% both_mRNA_positive_FC_gene, "positive", "negative")) %>%
  filter(Target_Entrez_Gene_ID %in% c(both_mRNA_negative_FC_gene, both_mRNA_positive_FC_gene)) %>%
  arrange(Target_Genome_Symbol) %>%
  write_tsv("gene_lists/miRNA_targets.tsv")