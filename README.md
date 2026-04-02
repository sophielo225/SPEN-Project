README
================
**About this project**

Solid pseudopapillary epithelial neoplasms (SPEN) are rare pancreatic tumors that predominantly affect young women and account for fewer than 3% of exocrine pancreatic neoplasms. Although mutations in the CTNNB1 gene and aberrant activation of the Wnt/β-catenin signaling pathway are known to play a role in tumor development, the downstream molecular mechanisms remain poorly understood. To address this gap, we performed an integrative transcriptomic analysis using publicly available datasets for SPEN. This includes two messenger RNA (mRNA) expression datasets and two microRNA (miRNA) expression datasets. While mRNAs reflect gene expression levels, miRNAs are small non-coding RNAs that regulate mRNA stability and translation, often playing a crucial role in tumor growth and suppression. These datasets were generated in different laboratories spanning 3 countries, using legacy experimental platforms. Thus, as a first step, we reprocessed the data using R/Bioconductor packages and recent annotations. We background corrected, normalized, and annotated the data before performing a differential expression analysis for each dataset. We then compared expression patterns across datasets to identify genes and miRNAs that were differentially expressed in multiple datasets with consistent directions of fold change, and performed gene-set analyses that incorporated known miRNA-mRNA interactions to identify coordinated regulatory signals associated with SPEN. This integrative approach should increase statistical power and enable the discovery of more consistent and biologically meaningful results that may be missed in single-dataset analyses. Our analysis highlights genes, miRNAs, molecular pathways, and regulatory interactions that may be involved in SPEN tumorigenesis and suggests mechanisms beyond CTNNB1 mutations. These findings provide a foundation for future biological studies and illustrate the value of reproducible data integration pipelines for extracting biological insights from heterogeneous public datasets.

**Information about getting the raw data**

We acquired datasets [GSE43795](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43795), [GSE43796](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43796), and [GSE140719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140719) directly from Gene Expression Omnibus (GEO). Because of file-formatting problems, we were unable to use the raw data directly for [E-MEXP-1914](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MEXP-1914). We processed it by manually downloading a binary R object posted on ArrayExpress, importing it into R, and then exporting it to the RDS format. That file is included here and is used in our analyses.

**Instructions on executing the code in this repository**

Please ensure that you install the latest version of [R](https://cran.rstudio.com/) and [RStudio](https://posit.co/download/rstudio-desktop/).
The code in this repository makes it possible to retrieve, process, and analyze the data from our analysis. To do so, please complete these steps:

1. Run `install.R` to install all the required packages
2. Run `SPEN_miRNA_analysis.R`. It contains code for processing the GSE43796 and GSE140719 datasets (miRNA).
3. Run `SPEN_mRNA_analysis.R`. It contains code for processing the GSE43795 and E-MEXP-1914 datasets (mRNA).
4. Run `SPEN_combined_analysis.R`. It contains code for doing the combined analysis for all four datasets.
5. Run `plot.R` to produce the graphs and tables.
