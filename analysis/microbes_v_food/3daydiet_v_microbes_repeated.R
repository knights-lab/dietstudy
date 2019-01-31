# plot food v microbe per person
setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_3day/")

# load datlist (out put from running correlations)
load(file = "../../../../data/food_v_microbes_per_person_3day.RData")

# write out a file showing all the significant correlations from each person
sigs <- lapply(datlist, function(x) subset(x, fdr_pval < 0.25))
allsigs <- do.call("rbind", sigs)

allsigs_names <- allsigs[colnames(allsigs) %in% c("Food", "Taxa")]

allsigs_names_pairs <- paste(allsigs_names$Food, allsigs_names$Taxa)

mytable<-as.data.frame(table(allsigs_names_pairs))

# look at the one that shows up 4 times
veg_lachno <- allsigs[grep("Other_vegetables_raw", allsigs$Food),]
veg_lachno <- veg_lachno[grep("f__Lachnospiraceae;NA;NA", veg_lachno$Taxa),]

ofinterest <- allsigs_names[duplicated(allsigs_names),]
ofinterest <- unique(ofinterest)




# Fix labeling
ofinterest$Food <- gsub(".*L3_", "", ofinterest$Food)
ofinterest$Food <- gsub("_", " ", ofinterest$Food)
ofinterest$Taxa <- gsub(".*;s__", "", ofinterest$Taxa)
ofinterest$Taxa <- gsub("?.*g__", "Uncl. Genus ", ofinterest$Taxa)
ofinterest$Taxa <- gsub("?.*f__", "Uncl. Family ", ofinterest$Taxa)
ofinterest$Taxa <- gsub("?.*o__", "Uncl. Order ", ofinterest$Taxa)
ofinterest$Taxa <- gsub("?.*p__", "Uncl. Phylum", ofinterest$Taxa)
ofinterest$Taxa <- gsub(";NA", "", ofinterest$Taxa)
ofinterest$Taxa <- gsub("_", " ", ofinterest$Taxa)

# change column names for supplemental table
colnames(ofinterest) <- c("Level 3 Food Groups", "Species or Higher Taxonomy")

write.table(ofinterest, file = "../../../../output/Table S1.txt", 
            sep = "\t", 
            col.names = T, 
            row.names = F,
            quote = F)

  meansigs <- NULL

for (i in 1:length(sigs)) {
  meansigs[i] <- dim(sigs[[i]])[1]
}

mean(meansigs)
sd(meansigs)

sum(meansigs>0)

