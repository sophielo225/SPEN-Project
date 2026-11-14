getmRNAExpressionPlotData <- function(exprs, labels, entrez_ids) {
  exprs = filter(exprs, Entrez_Gene_ID %in% entrez_ids) %>%
    pivot_longer(-Entrez_Gene_ID, names_to = "Sample_ID", values_to = "Value") %>%
    inner_join(labels) %>%
    inner_join(gene_annotations)

  return(exprs)
}

plotmRNAExpressionData <- function(data, title, out_file_path) {
  set.seed(0)
  
  ggplot(data, aes(x = Label, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Label), size = 3, alpha = 0.6) +
    theme_bw(base_size = 16) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(color = "") +
    facet_wrap(vars(HGNC_Symbol))
  
  ggsave(out_file_path, width = 6.5, height = 6)
}

# Read previously saved data
GSE43795_exprs = read_tsv("GSE43795_exprs.tsv.gz") %>%
  mutate(Entrez_Gene_ID = as.character(Entrez_Gene_ID))
GSE43795_labels = read_tsv("GSE43795_labels.tsv.gz")
e_mexp_1914_exprs = read_tsv("E-MEXP-1914_exprs.tsv.gz") %>%
  mutate(Entrez_Gene_ID = as.character(Entrez_Gene_ID))

plotmRNAExpressionData(getmRNAExpressionPlotData(GSE43795_exprs, GSE43795_labels, "4609"), "GSE43795 - underexpressed", "figures/GSE43795_mRNA_underexpressed.pdf")
plotmRNAExpressionData(getmRNAExpressionPlotData(GSE43795_exprs, GSE43795_labels, c("867", "9425", "23035", "5590", "9839")), "GSE43795 - overexpressed", "figures/GSE43795_mRNA_overexpressed.pdf")