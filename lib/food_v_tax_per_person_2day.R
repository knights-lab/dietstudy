#TODO: add a column for the person that the correlations come from for the big list of dataframes


# Make 2-day summary food table for the whole study

# read in individual summarized food tables for each person and combine them
setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_2day/")

# read all food files into a list
temp <- list.files(pattern = "*_food.txt")
foods_2d <- lapply(temp, function(x) read.delim(x, row.names = "taxonomy"))
foods_2d <- lapply(foods_2d, function(x) x[!names(x) == "X.FoodID"])
names(foods_2d) <- gsub("_food.txt", "", temp)

split = strsplit(rownames(foods_2d[[1]]),";")
foodStrings = sapply(split,function(x) paste(x[1:2],collapse=";"))
foods_2d <- lapply(foods_2d, function(x) rowsum(x,foodStrings))

# read in tax files into a list
temp <- list.files(pattern = "*_tax.txt")
tax_2d <- lapply(temp, function(x) read.delim(x, row.names = 1))
names(tax_2d) <- gsub("_tax.txt", "", temp)

# identify people with < 5 time points in the microbiome data
drops <- which(lapply(tax_2d, function(x) dim(x)[2])<=5)

# drop the people with too few time points
foods_2d <- foods_2d[!names(foods_2d) %in% names(drops)]
tax_2d <- tax_2d[!names(tax_2d) %in% names(drops)]



########### function to run correlations for each person ##############

run_correlations <- function(food.list, tax.list) {
  
  for (z in 1:32){
    
    x <- t(food.list[[z]])
    y <- t(tax.list[[z]])
   
    dat <- NULL
    for(i in colnames(x)){ # food
      for(j in colnames(y)){ # tax
        a<-as.numeric(x[,i])
        b<-as.numeric(y[,j])
        tmp <- c(i,
                 j,
                 suppressWarnings(cor(a,b, method = "spearman", use = "everything")), 
                 suppressWarnings(cor.test(a,b,method = "spearman")$p.value),
                 z)
        if(is.null(dat)){
          dat <- tmp
          }
        else{
          dat<-rbind(dat,tmp)
        }
      }
    }
    dat<- data.frame(row.names = NULL, dat, stringsAsFactors = FALSE)
    colnames(dat)<-c("Food","Taxa","Correlation","Pvalue", "ID")
    dat$Correlation <- as.numeric(dat$Correlation)
    dat$Pvalue <- as.numeric(dat$Pvalue)
    dat$fdr_pval <- (p.adjust(dat$Pvalue, 'fdr'))
    dat$fdr_pval <- ifelse(dat$Correlation == 1, 1, dat$fdr_pval)
    dat$Significance<-cut(dat$fdr_pval, breaks=c(-Inf, 0.25, Inf), label=c( "*", ""))
    
    # save this list of all values for future
    datlist[[z]] <<- dat
    
    # # prune for complete cases with significant correlations
    # dat.c<-dat[complete.cases(dat),]
    # taxkeep<-unique(dat.c$Taxa[dat.c$fdr_pval<=0.25])
    # 
    # taxassign<-taxkeep
    # dat.w<-dat.c[dat.c$Taxa %in% taxassign,] 
    # 
    # foodkeep <- unique(dat.w$Food[dat.w$fdr_pval<=0.25])
    # dat.m<-dat.w[dat.w$Food %in% foodkeep,]
    # 
    # 
    # # Fix labeling
    # dat.m$Food <- gsub(".*L3_", "", dat.m$Food)
    # dat.m$Food <- gsub("_", " ", dat.m$Food)
    # dat.m$Taxa <- gsub(".*;s__", "", dat.m$Taxa)
    # dat.m$Taxa <- gsub("?.*g__", "Uncl. Genus ", dat.m$Taxa)
    # dat.m$Taxa <- gsub("?.*f__", "Uncl. Family ", dat.m$Taxa)
    # dat.m$Taxa <- gsub("?.*o__", "Uncl. Order ", dat.m$Taxa)
    # dat.m$Taxa <- gsub("?.*p__", "Uncl. Phylum", dat.m$Taxa)
    # dat.m$Taxa <- gsub(";NA", "", dat.m$Taxa)
    # dat.m$Taxa <- gsub("_", " ", dat.m$Taxa)
    # 
    # 
    # # #plot the correlations for the collapsed levels
    # print(ggplot(data = dat.m, aes(x=Food, y=Taxa, fill=Correlation)) +
    #         geom_tile(color = "white") +
    #         # scale_fill_gradient2(low="#5e3c99", mid="#f7f7f7", high="#e66101", midpoint = 0, limit = c(-0.95,0.95), space = "Lab", name = "Spearman") +
    #         scale_fill_gradient2(low="purple", mid="white", high="green", midpoint = 0, limit = c(-.9999,.9999), space = "Lab", name = "Spearman") +
    #         geom_text(data = dat.m, aes(label=Significance), color="black", size=2) +
    #         #theme_minimal() +
    #         theme(axis.text.x = element_text(angle = 45, vjust = 1,
    #                                          size = 7, hjust = 1),
    #               axis.text.y = element_text(size = 7),
    #               axis.title = element_blank()) +
    #         #labs(x = "Species", y = "Species") +
    #         theme(strip.text.y = element_text(angle = 0, face = "italic"),
    #               strip.background = element_rect(color="grey", fill = "white")) +
    #         coord_fixed() +
    #         ggtitle(z))
    
  }
}

