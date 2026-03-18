# Install required packages if they are not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# CRAN packages
cran_packages <- c("tidyverse", "R.utils")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

# Bioconductor packages
bioc_packages <- c("GEOquery", "limma", "illuminaHumanv4.db", "hgug4110b.db", 
                   "oligo", "pd.mirna.4.0", "multiMiR", "miRBaseConverter")
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE)
    }
}
