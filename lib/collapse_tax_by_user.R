homedir <- "/Users/abby/Documents/Projects/dietstudy/"

require(robCompositions)
require(tibble)

setwd(dir = homedir)

# load map
tax_map <- read.table("data/maps/taxonomy_norm_map.txt", sep = "\t", header = T, comment = "")

# This file is pre-processed in clean taxonomy
tax <- read.delim("data/processed_tax/taxonomy_preprocessed.txt", row = 1, stringsAsFactors = F)

# Now collapse by UserName (UN)
UN_tax <- data.frame(matrix(nrow=nrow(tax), ncol=0))
rownames(UN_tax) <- rownames(tax)

# collapse by username (sum raw counts across the rows)
for (i in unique(tax_map$UserName)) {
  submap <- tax_map[tax_map$UserName == i,]
  subtax <- tax[,colnames(tax) %in% submap$X.SampleID]
  mytax <- as.data.frame(rowSums(subtax))
  colnames(mytax) <- i
  UN_tax <- cbind(UN_tax, mytax)
}

# Remove vary rare taxa from the tax table
#normalize
UN_tax_norm <- sweep(UN_tax, 2, colSums(UN_tax), "/")
#drop low abundance
UN_tax_norm_dla <- UN_tax_norm[rowMeans(UN_tax_norm) >= 0.0001,] # keeps 290 species-level anotations

# now limit tax to this
UN_tax <- UN_tax[rownames(UN_tax) %in% rownames(UN_tax_norm_dla),]

# Summarizing at different levels
split <- strsplit(rownames(UN_tax),";")             # Split and rejoin on lv7 to get species level

# Species
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
UN_tax_s <- rowsum(UN_tax,taxaStrings)            
rownames(UN_tax_s) = sapply(strsplit(rownames(UN_tax_s),";"),function(x) paste(x[1:7],collapse=";"));

# Genus
taxaStrings <- sapply(split,function(x) paste(x[1:6],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
UN_tax_g <- rowsum(UN_tax,taxaStrings)    
rownames(UN_tax_g) = sapply(strsplit(rownames(UN_tax_g),";"),function(x) paste(x[1:6],collapse=";"));

# Family
taxaStrings <- sapply(split,function(x) paste(x[1:5],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
UN_tax_f <- rowsum(UN_tax,taxaStrings) 
rownames(UN_tax_f) = sapply(strsplit(rownames(UN_tax_f),";"),function(x) paste(x[1:5],collapse=";"));


# Order
taxaStrings <- sapply(split,function(x) paste(x[1:4],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
UN_tax_o <- rowsum(UN_tax,taxaStrings) 
rownames(UN_tax_o) = sapply(strsplit(rownames(UN_tax_o),";"),function(x) paste(x[1:4],collapse=";"));

# Class
taxaStrings <- sapply(split,function(x) paste(x[1:3],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
UN_tax_c <- rowsum(UN_tax,taxaStrings) 
rownames(UN_tax_c) = sapply(strsplit(rownames(UN_tax_c),";"),function(x) paste(x[1:3],collapse=";"));

# Phylum
taxaStrings <- sapply(split,function(x) paste(x[1:2],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
UN_tax_p <- rowsum(UN_tax,taxaStrings) 
rownames(UN_tax_p) = sapply(strsplit(rownames(UN_tax_p),";"),function(x) paste(x[1:2],collapse=";"));



# create a normalized taxa tables at each level to use with some downstream analysis
##### Make relative abundance table for plotting bar graphs ####
# Normalize
UN_tax_norm_s <- sweep(UN_tax_s, 2, colSums(UN_tax_s), "/")
UN_tax_norm_g <- sweep(UN_tax_g, 2, colSums(UN_tax_g), "/")
UN_tax_norm_f <- sweep(UN_tax_f, 2, colSums(UN_tax_f), "/")
UN_tax_norm_o <- sweep(UN_tax_o, 2, colSums(UN_tax_o), "/")
UN_tax_norm_c <- sweep(UN_tax_c, 2, colSums(UN_tax_c), "/")
UN_tax_norm_p <- sweep(UN_tax_p, 2, colSums(UN_tax_p), "/")

# Sort by average abundance
UN_tax_norm_s <- UN_tax_norm_s[order(rowMeans(UN_tax_norm_s),decreasing=T),]
UN_tax_norm_g <- UN_tax_norm_g[order(rowMeans(UN_tax_norm_g),decreasing=T),]
UN_tax_norm_f <- UN_tax_norm_f[order(rowMeans(UN_tax_norm_f),decreasing=T),]
UN_tax_norm_o <- UN_tax_norm_o[order(rowMeans(UN_tax_norm_o),decreasing=T),]
UN_tax_norm_c <- UN_tax_norm_c[order(rowMeans(UN_tax_norm_c),decreasing=T),]
UN_tax_norm_p <- UN_tax_norm_p[order(rowMeans(UN_tax_norm_p),decreasing=T),]

# Drop low abundance taxa at each level
UN_tax_norm_s <- UN_tax_norm_s[rowMeans(UN_tax_norm_s) >= 0.0001,] # gets a little smaller
UN_tax_norm_g <- UN_tax_norm_g[rowMeans(UN_tax_norm_g) >= 0.0001,] # gets a little smaller
UN_tax_norm_f <- UN_tax_norm_f[rowMeans(UN_tax_norm_f) >= 0.0001,] # gets a little smaller
UN_tax_norm_o <- UN_tax_norm_o[rowMeans(UN_tax_norm_o) >= 0.0001,] # gets a little smaller
UN_tax_norm_c <- UN_tax_norm_c[rowMeans(UN_tax_norm_c) >= 0.0001,] # gets a little smaller
UN_tax_norm_p <- UN_tax_norm_p[rowMeans(UN_tax_norm_p) >= 0.0001,] # doesn't change

### Make clr-adjusted taxa table
# # 1. start with table of counts without bad samples or low abundance taxa
# # to limit to just the taxa that remain the the tax_norm table
# UN_tax_clr_s <- UN_tax_s[rownames(UN_tax_s) %in% rownames(UN_tax_norm_s),]
# UN_tax_clr_g <- UN_tax_g[rownames(UN_tax_g) %in% rownames(UN_tax_norm_g),]
# UN_tax_clr_f <- UN_tax_f[rownames(UN_tax_f) %in% rownames(UN_tax_norm_f),]
# UN_tax_clr_o <- UN_tax_o[rownames(UN_tax_o) %in% rownames(UN_tax_norm_o),]
# UN_tax_clr_c <- UN_tax_c[rownames(UN_tax_c) %in% rownames(UN_tax_norm_c),]
# UN_tax_clr_p <- UN_tax_p[rownames(UN_tax_p) %in% rownames(UN_tax_norm_p),]

# 1. start with table of relative abundances
UN_tax_clr_s <- UN_tax_norm_s
UN_tax_clr_g <- UN_tax_norm_g
UN_tax_clr_f <- UN_tax_norm_f
UN_tax_clr_o <- UN_tax_norm_o
UN_tax_clr_c <- UN_tax_norm_c
UN_tax_clr_p <- UN_tax_norm_p


# 2. interpolate 0s
mynames <- colnames(UN_tax)    # store colnames for later

# # change int > num
# UN_tax_clr_s <- UN_tax_clr_s*1.0      
# UN_tax_clr_g <- UN_tax_clr_g*1.0  
# UN_tax_clr_f <- UN_tax_clr_f*1.0  
# UN_tax_clr_o <- UN_tax_clr_o*1.0  
# UN_tax_clr_c <- UN_tax_clr_c*1.0  
# UN_tax_clr_p <- UN_tax_clr_p*1.0 

# transpose for imputation
transposed_s <- t(UN_tax_clr_s)
transposed_g <- t(UN_tax_clr_g)
transposed_f <- t(UN_tax_clr_f)
transposed_o <- t(UN_tax_clr_o)
transposed_c <- t(UN_tax_clr_c)
transposed_p <- t(UN_tax_clr_p)

# #this is the imputation step and it takes ages to finish, so just run when needed
myimpR_s = impRZilr(transposed_s, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_s)), verbose = T)
myimpR_g = impRZilr(transposed_g, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_g)), verbose = T)
myimpR_f = impRZilr(transposed_f, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_f)), verbose = T)
myimpR_o = impRZilr(transposed_o, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_o)), verbose = T)
myimpR_c = impRZilr(transposed_c, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_c)), verbose = T)
myimpR_p = impRZilr(transposed_p, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_p)), verbose = T)

# save as rdata to use again later without re-running imputation
save(myimpR_s, file = "data/processed_UN_tax/myimpR_UN_s.rdata")
save(myimpR_g, file = "data/processed_UN_tax/myimpR_UN_g.rdata")
save(myimpR_f, file = "data/processed_UN_tax/myimpR_UN_f.rdata")
save(myimpR_o, file = "data/processed_UN_tax/myimpR_UN_o.rdata")
save(myimpR_c, file = "data/processed_UN_tax/myimpR_UN_c.rdata")
save(myimpR_p, file = "data/processed_UN_tax/myimpR_UN_p.rdata")

# load files for use
load(file = "data/processed_UN_tax/myimpR_UN_s.rdata")
load(file = "data/processed_UN_tax/myimpR_UN_g.rdata")
load(file = "data/processed_UN_tax/myimpR_UN_f.rdata")
load(file = "data/processed_UN_tax/myimpR_UN_o.rdata")
load(file = "data/processed_UN_tax/myimpR_UN_c.rdata")
load(file = "data/processed_UN_tax/myimpR_UN_p.rdata")

# 3. clr tranformation
UN_tax_clr_s = t(cenLR(myimpR_s$x)$x.clr)  # clr tranformation of imputed table
UN_tax_clr_g = t(cenLR(myimpR_g$x)$x.clr)
UN_tax_clr_f = t(cenLR(myimpR_f$x)$x.clr)
UN_tax_clr_o = t(cenLR(myimpR_o$x)$x.clr)
UN_tax_clr_c = t(cenLR(myimpR_c$x)$x.clr)
UN_tax_clr_p = t(cenLR(myimpR_p$x)$x.clr)

# re-add names
colnames(UN_tax_clr_s) <- mynames
colnames(UN_tax_clr_g) <- mynames
colnames(UN_tax_clr_f) <- mynames
colnames(UN_tax_clr_o) <- mynames
colnames(UN_tax_clr_c) <- mynames
colnames(UN_tax_clr_p) <- mynames

# 4. now good to go for export/use with stats
UN_tax_clr_s <- as.data.frame(UN_tax_clr_s)
UN_tax_clr_g <- as.data.frame(UN_tax_clr_g)
UN_tax_clr_f <- as.data.frame(UN_tax_clr_f)
UN_tax_clr_o <- as.data.frame(UN_tax_clr_o)
UN_tax_clr_c <- as.data.frame(UN_tax_clr_c)
UN_tax_clr_p <- as.data.frame(UN_tax_clr_p)

#### Make counts table for alpha diversity from normalized table #####
# multiply re-normalized reduced table by a factor to remove decimals and round
UN_tax_counts_s <- round(sweep(UN_tax_norm_s, 2, colSums(UN_tax_norm_s),'/')*10000)
UN_tax_counts_g <- round(sweep(UN_tax_norm_g, 2, colSums(UN_tax_norm_s),'/')*10000)
UN_tax_counts_f <- round(sweep(UN_tax_norm_f, 2, colSums(UN_tax_norm_s),'/')*10000)
UN_tax_counts_o <- round(sweep(UN_tax_norm_o, 2, colSums(UN_tax_norm_s),'/')*10000)
UN_tax_counts_c <- round(sweep(UN_tax_norm_c, 2, colSums(UN_tax_norm_s),'/')*10000)
UN_tax_counts_p <- round(sweep(UN_tax_norm_p, 2, colSums(UN_tax_norm_s),'/')*10000)

# prep for export
UN_tax_norm_s <- rownames_to_column(UN_tax_norm_s, var = "#taxonomy")
UN_tax_norm_g <- rownames_to_column(UN_tax_norm_g, var = "#taxonomy")
UN_tax_norm_f <- rownames_to_column(UN_tax_norm_f, var = "#taxonomy")
UN_tax_norm_o <- rownames_to_column(UN_tax_norm_o, var = "#taxonomy")
UN_tax_norm_c <- rownames_to_column(UN_tax_norm_c, var = "#taxonomy")
UN_tax_norm_p <- rownames_to_column(UN_tax_norm_p, var = "#taxonomy")

UN_tax_counts_s <- rownames_to_column(UN_tax_counts_s, var = "#taxonomy")
UN_tax_counts_g <- rownames_to_column(UN_tax_counts_g, var = "#taxonomy")
UN_tax_counts_f <- rownames_to_column(UN_tax_counts_f, var = "#taxonomy")
UN_tax_counts_o <- rownames_to_column(UN_tax_counts_o, var = "#taxonomy")
UN_tax_counts_c <- rownames_to_column(UN_tax_counts_c, var = "#taxonomy")
UN_tax_counts_p <- rownames_to_column(UN_tax_counts_p, var = "#taxonomy")


# prep clr tables for export
UN_tax_clr_s <- rownames_to_column(UN_tax_clr_s, var = "#taxonomy")
UN_tax_clr_g <- rownames_to_column(UN_tax_clr_g, var = "#taxonomy")
UN_tax_clr_f <- rownames_to_column(UN_tax_clr_f, var = "#taxonomy")
UN_tax_clr_o <- rownames_to_column(UN_tax_clr_o, var = "#taxonomy")
UN_tax_clr_c <- rownames_to_column(UN_tax_clr_c, var = "#taxonomy")
UN_tax_clr_p <- rownames_to_column(UN_tax_clr_p, var = "#taxonomy")

# export normalized tables
write.table(UN_tax_norm_s, file="data/processed_UN_tax/UN_taxonomy_norm_s.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_norm_g, file="data/processed_UN_tax/UN_taxonomy_norm_g.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_norm_f, file="data/processed_UN_tax/UN_taxonomy_norm_f.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_norm_o, file="data/processed_UN_tax/UN_taxonomy_norm_o.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_norm_c, file="data/processed_UN_tax/UN_taxonomy_norm_c.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_norm_p, file="data/processed_UN_tax/UN_taxonomy_norm_p.txt", quote=F, sep="\t", row.names = F, col.names = T)

# export counts tables (really probably only need counts at the species level)
write.table(UN_tax_counts_s, file = "data/processed_UN_tax/UN_taxonomy_counts_s.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(UN_tax_counts_g, file = "data/processed_UN_tax/UN_taxonomy_counts_g.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(UN_tax_counts_f, file = "data/processed_UN_tax/UN_taxonomy_counts_f.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(UN_tax_counts_o, file = "data/processed_UN_tax/UN_taxonomy_counts_o.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(UN_tax_counts_c, file = "data/processed_UN_tax/UN_taxonomy_counts_c.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(UN_tax_counts_p, file = "data/processed_UN_tax/UN_taxonomy_counts_p.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# export clr transformed taxonomy tables for each level
write.table(UN_tax_clr_s, file = "data/processed_UN_tax/UN_taxonomy_clr_s.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_clr_g, file = "data/processed_UN_tax/UN_taxonomy_clr_g.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_clr_f, file = "data/processed_UN_tax/UN_taxonomy_clr_f.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_clr_o, file = "data/processed_UN_tax/UN_taxonomy_clr_o.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_clr_c, file = "data/processed_UN_tax/UN_taxonomy_clr_c.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(UN_tax_clr_p, file = "data/processed_UN_tax/UN_taxonomy_clr_p.txt", quote = F, sep="\t", row.names = F, col.names = T)

