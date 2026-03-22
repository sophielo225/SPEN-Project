README
================
**About this project**

Solid pseudopapillary epithelial neoplasms (SPEN) are rare pancreatic tumors that predominantly affect young women and account for fewer than 3% of exocrine pancreatic neoplasms. Although mutations in the CTNNB1 gene and aberrant activation of the Wnt/β-catenin signaling pathway are known to play a role in tumor development, the downstream molecular mechanisms remain poorly understood. To address this gap, we performed an integrative transcriptomic analysis using publicly available datasets for this disease: two messenger RNA (mRNA) expression datasets and two microRNA (miRNA) expression datasets. While mRNAs reflect gene expression levels, miRNAs are small non-coding RNAs that regulate mRNA stability and translation, often playing a crucial role in tumor growth and suppression. These datasets were generated in different laboratories spanning 3 countries, using legacy experimental platforms. Thus, as a first step, we reprocessed the data using R/Bioconductor packages. We background corrected, normalized, and annotated the data before performing a differential expression analysis for each dataset. We then compared expression patterns across datasets to identify genes and miRNAs that were differentially expressed in multiple datasets with consistent directions of fold change, and performed gene-set analyses that incorporated known miRNA-mRNA interactions to identify coordinated regulatory signals associated with SPEN. This integrative approach increases statistical power and enables the discovery of more consistent and biologically meaningful results that may be missed in single-dataset analyses. Our analysis highlights genes, miRNAs, molecular pathways and regulatory interactions that may be involved in SPEN tumorigenesis and suggests mechanisms beyond CTNNB1 mutations. These findings provide a foundation for future biological studies and illustrate the value of reproducible data integration pipelines for extracting biological insights from heterogeneous public datasets. <br><br>

**Information about getting the raw data**

GSE43795, GSE43796, and GSE140719 were acquired from Gene Expression Omnibus (GEO). Because of versioning problem, E-MEXP-1914 was accessed and processed by donwloading its RDS object on ArrayExpress. All GEO datasets are able to be directly downloaded through the codes, while E-MEXP-1914 RDS object is included in this repo and used for further processing. <br><br>

**Instructions on running files in this repo**

To be able to regenerate the results we got in this study, follow the steps below:
1. Run install.R file to install all the required packages
2. Run SPEN_miRNA_analysis.R file. It contains codes for processing GSE43796 and GSE140719 datasets.
3. Run SPEN_mRNA_analysis.R file. It contains codes for processing GSE43795 and E-MEXP-1914 datasets.
4. Run SPEN_combined_analysis.R file. It contains codes for doing combined analysis for all four datasets.
5. Run plotting.R file to get graphs and tables.
