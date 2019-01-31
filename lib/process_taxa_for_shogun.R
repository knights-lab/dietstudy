library(tibble)

taxa <- read.delim(file = "~/Documents/Projects/dietstudy/data/shogun_results_dietstudy_170822/embalmer_taxatable.txt", comment = "", row.names = 1)

#### Manicure the samplenames, grab latest (mct study only!)
colnames(taxa) <- gsub(".S[0-9]+.R1.001","",colnames(taxa))      # Clean old plate IDs
taxa <- taxa[,order(colnames(taxa))]              # Sort nicely by sample ID
taxa <- taxa[,-(grep("L0",colnames(taxa))-1)]    # Keep new runs only
colnames(taxa) <- gsub(".S[0-9]+.L001.R1.001","",colnames(taxa)) # Clean new plate IDs

#### Drop blanks 
taxa <- taxa[,-(grep("Blank", colnames(taxa)))]


### Drop low-read samples
taxa <- taxa[, colSums(taxa) > 20000]

taxa <- rownames_to_column(taxa, var = "#OTU")

write.table(taxa, file = "~/Documents/Projects/dietstudy/data/shogun_results_cleaned/embalmer_taxatable_clean.txt", 
            sep = "\t", 
            quote = F,
            col.names = T,
            row.names = F)
