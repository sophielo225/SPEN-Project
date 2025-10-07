library(tidyverse)

# Read files without the headers and descriptions, only keep the expression data
GSE140719 = read_tsv("GSE140719_series_matrix.txt", skip = 71)

# Read reference file for GSE140719, also exclude the headers
GSE140719_ref = read_tsv("ID reference for GSE140719.txt", skip = 17) %>%
    rename(ID_REF = ID)

# Map reference data to GSE140719
GSE140719_annotated = left_join(GSE140719, GSE140719_ref, by = "ID_REF") %>%
    rename(sequence_type = 'Sequence Type')

# Filter out all non-miRNA entries
GSE140719_MIMAT_annotated = filter(GSE140719_annotated, sequence_type == "miRNA")

# Get MIMAT accession ID for GSE140719
GSE140719_MIMAT_ID = pull(GSE140719_MIMAT_annotated, Accession)


# Read GSE43796 file without the headers and descriptions, only keep the expression data
GSE43796 = head(read_tsv("GSE43796_series_matrix.txt", skip = 58), -1)

# Read reference file for GSE43796
GSE43796_ref = read_tsv("GSE43796_ID_Ref.txt")

# The MIMAT accession ID is the stable, unique ID across miRBase system
GSE43796_ref_MIMAT = read_tsv("GSE43796_ID_Ref.txt") %>%
    select(ID, miRNA_ID, ACCESSION_STRING) %>%
    rename(ID_REF = ID) %>%
    mutate(MIMAT = str_extract(ACCESSION_STRING, "MIMAT\\d+"))

# Map reference MIMAT accession ID to GSE43796
GSE43796_annotated = left_join(GSE43796, GSE43796_ref_MIMAT, by = "ID_REF")

# Get MIMAT accession ID for GSE43796
GSE43796_MIMAT_ID = pull(GSE43796_annotated, MIMAT)

# Get MIMAT accession ID that are both in GSE43796 and GSE140719
# There are 1197 mutual miRNA IDs
mutual_MIMAT_ID = intersect(GSE140719_MIMAT_ID, GSE43796_MIMAT_ID)

# Filter the data sets to only contain the mutual MIMAT ID
GSE140719_filtered = filter(GSE140719_MIMAT_annotated, Accession %in% mutual_MIMAT_ID) %>%
    select(Accession, GSM4182441:GSM4182451) %>%
    rename(MIMAT_ID = Accession)

GSE43796_filtered = filter(GSE43796_annotated, MIMAT %in% mutual_MIMAT_ID) %>%
    select(MIMAT, GSM1071387:GSM1071417) %>%
    rename(MIMAT_ID = MIMAT)

# Since GSE43796 has duplicate expression level for some miRNA ID, 
# we need to pivot the tibble and calculate the average level if the miRNA ID has multiple expression levels
GSE43796_pivoted = pivot_longer(GSE43796_filtered, GSM1071387:GSM1071417, names_to = "Sample_name", values_to = "Expression_level") %>%
    group_by(MIMAT_ID, Sample_name) %>%
    summarize(mean_expression = mean(Expression_level)) %>%
    pivot_wider(names_from = "Sample_name", values_from = "mean_expression")

# Filter the data sets to contain only the samples of interest
# For GSE43796, GSM1071387 to GSM1071400 are SPN tumor samples, and GSM1071413 to GSM1071417 are normal samples
# For GSE140719, localized tumor samples are: GSM4182441, GSM4182443, GSM4182444, GSM4182447, GSM4182451
# and normal samples are: GSM4182442, GSM4182446, GSM4182448
GSE43796_SPN = select(GSE43796_pivoted, GSM1071387:GSM1071400, GSM1071413:GSM1071417)
GSE140719_SPN = select(GSE140719_filtered, MIMAT_ID, GSM4182441, GSM4182443, GSM4182444, GSM4182447, 
                       GSM4182451, GSM4182442, GSM4182446, GSM4182448)

# Prepare data sets to perform Mann-Whitney test
# GSE43796
GSE43796_SPN = pivot_longer(
    GSE43796_SPN, cols = starts_with("GSM"),
    names_to = "sample_id",
    values_to = "expression"
)
sample_groups <- tibble(
    sample_id = c("GSM1071387","GSM1071388","GSM1071389","GSM1071390","GSM1071391",
                  "GSM1071392", "GSM1071393", "GSM1071394", "GSM1071395", "GSM1071396",
                  "GSM1071397", "GSM1071398", "GSM1071399", "GSM1071400", "GSM1071413",
                  "GSM1071414", "GSM1071415", "GSM1071416", "GSM1071417"),
    group = c(rep("tumor", 14), rep("normal", 5))
)

GSE43796_tidy <- GSE43796_SPN %>%
    left_join(sample_groups, by = "sample_id")

# GSE140719
GSE140719_SPN = pivot_longer(
    GSE140719_SPN, cols = starts_with("GSM"),
    names_to = "sample_id",
    values_to = "expression"
)
sample_groups <- tibble(
    sample_id = c("GSM4182441", "GSM4182443", "GSM4182444", "GSM4182447", "GSM4182451",
                  "GSM4182442", "GSM4182446", "GSM4182448"),
    group = c(rep("tumor", 5), rep("normal", 3))
)

GSE140719_tidy <- GSE140719_SPN %>%
    left_join(sample_groups, by = "sample_id")


# Perform Mann-Whitney Test & find q-value
GSE43796_result = GSE43796_tidy %>%
    group_by(MIMAT_ID) %>%
    summarise(p_value = wilcox.test(expression ~ group)$p.value) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))

GSE140719_result = GSE140719_tidy %>%
    group_by(MIMAT_ID) %>%
    summarise(p_value = wilcox.test(expression ~ group)$p.value) %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))

# Filter out statistically significant IDs that are in both data sets
# Set the threshold into 0.05
GSE43796_result %>%
    filter(q_value < 0.05) %>%
    pull(MIMAT_ID) -> GSE43796_MIMAT_sig
print(GSE43796_MIMAT_sig)

GSE140719_result %>%
    filter(q_value < 0.05) %>%
    pull(MIMAT_ID) -> GSE140719_MIMAT_sig
print(GSE140719_MIMAT_sig) # This is zero!! (why????)

mutual_MIMAT_sig = intersect(GSE43796_MIMAT_sig, GSE140719_MIMAT_sig)
print(mutual_MIMAT_sig)






