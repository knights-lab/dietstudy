homedir <- "/Users/abby/Documents/Projects/dietstudy/"
setwd(dir = homedir)

# Load tables from shogun
kegg <- read.delim("data/shogun_results_cleaned/shogun_function/filtered_functions/embalmer_taxatable_clean.species.kegg.txt", row = 1, stringsAsFactors = F)
modules <- read.delim("data/shogun_results_cleaned/shogun_function/filtered_functions/embalmer_taxatable_clean.species.kegg.modules.txt", row = 1, stringsAsFactors = F)
coverage <- read.delim("data/shogun_results_cleaned/shogun_function/filtered_functions/embalmer_taxatable_clean.species.kegg.modules.coverage.txt", row = 1, stringsAsFactors = F)

# make modules equal to modules * coverage
# check order
rownames(modules) == rownames(coverage)
colnames(modules) == colnames(coverage)
modules <- modules * coverage

rownames(modules) <- gsub(" ", "", rownames(modules))

# load map
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# Drop dropouts 
dropouts =  as.character(map[map$UserName %in% c("MCTs02", "MCTs17", "MCTs30"),"X.SampleID"])
kegg <- kegg[, !(colnames(kegg) %in% dropouts)]
modules <- modules[, !(colnames(modules) %in% dropouts)]

# Make shogun map
shogun_map <- map[map$X.SampleID %in% colnames(kegg),]


# limit modules to pathway
# import kegg module map
modules_map <- read.table("raw/functional/Kegg_ID_Map_plus.txt", header = T, sep = "\t", comment = "", stringsAsFactors = F)
path_modules_map <- modules_map[grep("Pathway module", modules_map$Reference_hierarchy),]

#limit to just those modules
pathway_modules <- modules[rownames(modules) %in% path_modules_map$Shogun_ID,]

# write this version to a file for use in other scripts
write.table(kegg, file="data/shogun_results_cleaned/shogun_no_dropouts/kegg.txt", quote=F, sep="\t", row.names = T, col.names = T)
write.table(modules, file="data/shogun_results_cleaned/shogun_no_dropouts/modules.txt", quote=F, sep="\t", row.names = T, col.names = T)
write.table(pathway_modules, file="data/shogun_results_cleaned/shogun_no_dropouts/pathway_modules.txt", quote=F, sep="\t", row.names = T, col.names = T)


# Summarize by UserName
kegg_smry <- data.frame(matrix(nrow=nrow(kegg), ncol=0))
modules_smry <- data.frame(matrix(nrow=nrow(modules), ncol = 0))

rownames(kegg_smry) <- rownames(kegg)
rownames(modules_smry) <- rownames(modules)

for (i in unique(shogun_map$UserName)){
  subgroup <- shogun_map[shogun_map$UserName == i,]
  # kegg
  keggsub <- kegg[,colnames(kegg) %in% subgroup$X.SampleID]
  tmp <- as.data.frame(rowMeans(keggsub))                    # need to choose if means or sums is better for analysis
  colnames(tmp) <- i
  kegg_smry <- cbind(kegg_smry,tmp)
  # modules
  modulessub <- modules[,colnames(modules) %in% subgroup$X.SampleID]
  tmp1 <- as.data.frame(rowMeans(modulessub))
  colnames(tmp1) <- i
  modules_smry <- cbind(modules_smry, tmp1)
}

#limits modules_smry to just pathway modules
path_mod_smry <- modules_smry[rownames(modules_smry) %in% path_modules_map$Shogun_ID,]


# write username summaries to files
write.table(kegg_smry, file="data/shogun_results_cleaned/shogun_no_dropouts/kegg_smry.txt", quote=F, sep="\t", row.names = T, col.names = T)
write.table(modules_smry, file="data/shogun_results_cleaned/shogun_no_dropouts/modules_smry.txt", quote=F, sep="\t", row.names = T, col.names = T)
write.table(path_mod_smry, file="data/shogun_results_cleaned/shogun_no_dropouts/path_modules_smry.txt", quote=F, sep="\t", row.names = T, col.names = T)
