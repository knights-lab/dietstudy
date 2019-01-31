
setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_2day/")

load(file = "../../../../data/food_v_microbes_per_person_2day.RData")

# write out a file showing all the significant correlations from each person
sigs <- lapply(datlist, function(x) subset(x, fdr_pval < 0.25))
allsigs <- do.call("rbind", sigs)

allsigs_names <- allsigs[colnames(allsigs) %in% c("Food", "Taxa")]

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

meansigs <- NULL

for (i in 1:length(sigs)) {
  meansigs[i] <- dim(sigs[[i]])[1]
}

mean(meansigs)
sd(meansigs)
