getmiRNAExpressionPlotData <- function(exprs, labels, ids) {
  inner_join(exprs, miRNA_names, by = c("MIMAT_ID" = "Accession")) %>%
    mutate(miRNA_ID = TargetName) %>%
    dplyr::select(-MIMAT_ID, -TargetName) %>%
    dplyr::select(miRNA_ID, everything()) %>%
    filter(miRNA_ID %in% ids) %>%
    pivot_longer(
      -miRNA_ID,
      names_to = "Sample_ID",
      values_to = "Value"
    ) %>%
    inner_join(labels, by = "Sample_ID") %>%
    return()
}

plotmiRNAExpressionData <- function(data, title, out_file_path) {
  set.seed(0)
  
  ggplot(data, aes(x = Label, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Label), size = 3, alpha = 0.6) +
    theme_bw(base_size = 14) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(color = "") +
    facet_wrap(vars(miRNA_ID))
  
  ggsave(out_file_path, width = 6.5, height = 6)
}

GSE43796_exprs = read_tsv("GSE43796_exprs.tsv.gz")
GSE43796_labels = read_tsv("GSE43796_labels.tsv.gz")
GSE140719_exprs = read_tsv("GSE140719_exprs.tsv.gz")
GSE140719_labels = read_tsv("GSE140719_labels.tsv.gz")

plotmiRNAExpressionData(getmiRNAExpressionPlotData(GSE43796_exprs, GSE43796_labels, negative_FC_miRNA), "GSE43796 - underexpressed", "figures/GSE43796_miRNA_underexpressed.pdf")
plotmiRNAExpressionData(getmiRNAExpressionPlotData(GSE43796_exprs, GSE43796_labels, positive_FC_miRNA), "GSE43796 - overexpressed", "figures/GSE43796_miRNA_overexpressed.pdf")
plotmiRNAExpressionData(getmiRNAExpressionPlotData(GSE140719_exprs, GSE140719_labels, negative_FC_miRNA), "GSE140719 - underexpressed", "figures/GSE140719_miRNA_underexpressed.pdf")
plotmiRNAExpressionData(getmiRNAExpressionPlotData(GSE140719_exprs, GSE140719_labels, positive_FC_miRNA), "GSE140719 - overexpressed", "figures/GSE140719_miRNA_overexpressed.pdf")