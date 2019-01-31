#TODO: run collapse by correlation first before processing

homedir <- "/Users/abby/Documents/Projects/dietstudy/"

require(robCompositions)
require(tibble)

setwd(dir = homedir)

# load map
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# newest taxonomy file from October 2
tax <- read.delim("raw/hoho.tax.txt", row = 1, stringsAsFactors = F)

# Drop samples with low read counts (for original table this was >=20000)
# Keep everything after the bad high-read blank
tax <- tax[, colSums(tax) >=23500]

# Drop dropouts 
dropouts =  as.character(map[map$UserName %in% c("MCTs02", "MCTs17", "MCTs30"),"X.SampleID"])
tax <- tax[, !(colnames(tax) %in% dropouts)]

# see how many non-bacterial hits
tax_non_bact <- tax[-grep("k__Bacteria", rownames(tax)),]

# only a few, but a couple have decent counts in some people, so they may survive filtering and I will keep them in.

# DONT DO ANYMORE: remove non-bacteria taxa
# tax <- tax[grep("k__Bacteria", rownames(tax)),]

# make a limited map that matches all samples left
tax_map <- map[map$X.SampleID %in% colnames(tax),]


# Drop taxa that don't appear in at least 25% of people in the study
bugs_per_person <- list()

for (i in unique(tax_map$UserName)) {
  submap <- tax_map[tax_map$UserName == i,]
  subset <- tax[, colnames(tax) %in% submap$X.SampleID]
  subset <- subset[rowSums(subset > 0) > ncol(subset)/4, ] # these are the bugs present in this person at least 1/4 days in the study
  mybugs <- rownames(subset)
  bugs_per_person[[i]] <- rownames(subset)
}

n25 <- round(length(bugs_per_person) * 0.25) # 10% of people is about 8
nonuni <- unlist(bugs_per_person)
counts <- table(nonuni)
counts <- counts[counts > n25]
keep <- as.vector(names(counts))

tax <- tax[rownames(tax) %in% keep,]  # limit to taxa we care about
tax <- tax[rowSums(tax) > 0,]         # double check all 0 sum taxa are gone


# Remove vary rare taxa from the tax table
  #normalize
tax_norm <- sweep(tax, 2, colSums(tax), "/")
  #drop low abundance (this is more strict than previous filtering for other tasks)
tax_norm_dla <- tax_norm[rowMeans(tax_norm) >= 0.0005,] # keeps 153 species-level anotations

# now limit tax to this
tax <- tax[rownames(tax) %in% rownames(tax_norm_dla),]

# Summarizing at different levels
split <- strsplit(rownames(tax),";")             # Split and rejoin on lv7 to get species level

# Species
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_s <- rowsum(tax,taxaStrings)
rownames(tax_s) = sapply(strsplit(rownames(tax_s),";"),function(x) paste(x[1:7],collapse=";"));

# now we have about 115 species-level assignments


# create a normalized taxa tables at each level to use with some downstream analysis
##### Make relative abundance table for plotting bar graphs ####
# Normalize
tax_norm_s <- sweep(tax_s, 2, colSums(tax_s), "/")

# Sort by average abundance
tax_norm_s <- tax_norm_s[order(rowMeans(tax_norm_s),decreasing=T),]


### Make clr-adjusted taxa table
# 1. start with table of counts without bad samples or low abundance taxa
# to limit to just the taxa that remain the the tax_norm table
# columns are samples
tax_clr_s <- tax_s[rownames(tax_s) %in% rownames(tax_norm_s),]


# 2. interpolate 0s
mynames <- colnames(tax)    # store colnames (sample names) for later

# change int > num
tax_clr_s <- tax_clr_s*1.0      

# transpose for imputation
# columns are taxa
transposed_s <- t(tax_clr_s)


#this is the imputation step and it takes ages to finish, so just run when needed
#myimpR_s = impRZilr(transposed_s, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_s)), verbose = T)

# save as rdata to use again later without re-running imputation
#save(myimpR_s, file = "data/processed_tax/myimpR_s_for_corr.rdata")

# load files for use
load(file = "data/processed_tax/myimpR_s_for_corr.rdata")

# 3. clr tranformation
tax_clr_s = t(cenLR(myimpR_s$x)$x.clr)  # clr tranformation of imputed table

# re-add names
colnames(tax_clr_s) <- mynames

# 4. now good to go for export/use with stats
tax_clr_s <- as.data.frame(tax_clr_s)




# TODO:make a correlation collapsed table
######### collapse by correlation ##########

#Collapse species that are highly correlated
my_table <- tax_clr_s

#Get Correlations for module table
Cor_outs <- cor(t(my_table))

#view corrleations (red are highly correlated)
image(Cor_outs)

#Collapse by correlation
#Select rows where the correlation is > 0.95 
Cor_outs2 <- abs(Cor_outs) > 0.95
new_mytable <- my_table

for(i in 1:nrow(Cor_outs2)){
  for(j in i:nrow(Cor_outs2)){
    if(i==j){
      next
    }
    if(is.na(Cor_outs2[j,i])){
      next
    }
    if(Cor_outs2[j,i]){
      new_mytable[j,] <- new_mytable[j,] + new_mytable[i,]
      new_mytable[i,] <- rep(0, ncol(new_mytable))
    }
  }
}

#Collapse
my_table.new <- new_mytable[rowSums(new_mytable) > 0, ]
#sanity check
#sum(taxa_table.new) == sum(taxa_table)

#New correlation
cor_out3 <- cor(t(my_table.new))

#print new heatmap 
image(cor_out3)

# limit to highly correlated pathways (126)
keep <- rownames(my_table.new)
tax_clr_s_corr_colapse <- tax_clr_s[keep,]

##################################

#### Make counts table for alpha diversity from normalized table #####
# multiply re-normalized reduced table by a factor to remove decimals and round
tax_counts_s <- round(sweep(tax_norm_s, 2, colSums(tax_norm_s),'/')*med_depth)


# prep for export
tax_norm_s <- rownames_to_column(tax_norm_s, var = "#taxonomy")

tax_counts_s <- rownames_to_column(tax_counts_s, var = "#taxonomy")


# fix names on taxonomy map
colnames(tax_map)[1] <- "#SampleID"

# prep clr tables for export
tax_clr_s <- rownames_to_column(tax_clr_s, var = "#taxonomy")
tax_clr_s_corr_colapse <- rownames_to_column(tax_clr_s_corr_colapse, var = "#taxonomy")

# export normalized tables
write.table(tax_norm_s, file="data/processed_tax/taxonomy_norm_s_for_corr.txt", quote=F, sep="\t", row.names = F, col.names = T)

# export counts tables (really probably only need counts at the species level)
write.table(tax_counts_s, file = "data/processed_tax/taxonomy_counts_s_for_corr.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# export clr transformed taxonomy tables for each level
write.table(tax_clr_s, file = "data/processed_tax/taxonomy_clr_s_for_corr.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(tax_clr_s_corr_colapse, file = "data/processed_tax/taxonomy_clr_s_for_corr_collapsed.txt", quote = F, sep = "\t", row.names = F, col.names = T)
