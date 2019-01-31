
setwd(dir = "/Users/abby/Documents/Projects/dietstudy/")
require(ape)
require(tidyverse)

# load map
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID", "UserName","StudyDayNo","Sample.Week.Day")]


# load diet 3 day beta div for all people
diet3daydist <- read.delim("data/processed_food/3dayfood_beta/unweighted_unifrac_master3daydiet.txt", row.names = 1)

# get diet pcs
diet3daypc <- as.data.frame(pcoa(diet3daydist)$vectors)

# load mb table for everyone
mb <- read.delim("data/processed_tax/taxonomy_clr_s.txt", row = 1)

# get mb pcs
mbdist <- dist(t(mb))
mbpc <- as.data.frame(pcoa(mbdist)$vectors)

# limit mb and mbpc to just the people with diet data
mb <- mb[,colnames(mb) %in% rownames(diet3daypc)]
mbpc <- mbpc[rownames(mbpc) %in% rownames(diet3daypc),]

# limit diet to just the days with mb data
diet3daypc <- diet3daypc[rownames(diet3daypc) %in% rownames(mbpc),]

# limit map to just these people too
map <- map[map$X.SampleID %in% rownames(mbpc),]
map <- droplevels(map)

# limit both diet and mb pcs to the top 5 pcs
mbpc <- mbpc[,1:5]
diet3daypc <- diet3daypc[,1:5]

# fix naming
colnames(mbpc) <- paste0("Mb.Axis.", 1:5)
colnames(diet3daypc) <- paste0("Dt.Axis.", 1:5)

# clean up naming in the mb df
# fix the naming of each species
fix.names <- function(x) {
  rownames(x) <- gsub("?.*s__", "", rownames(x))
  rownames(x) <- gsub("?.*g__", "Uncl. Genus ", rownames(x))
  rownames(x) <- gsub("?.*f__", "Uncl. Family ", rownames(x))
  rownames(x) <- gsub("?.*o__", "Uncl. Order ", rownames(x))
  rownames(x) <- gsub("?.*p__", "Uncl. Phylum ", rownames(x))
  rownames(x) <- gsub("?.*k__", "Uncl. Kingdom ", rownames(x))
  rownames(x) <- gsub(";NA", "", rownames(x))
  rownames(x) <- gsub("\\[", "", rownames(x))
  rownames(x) <- gsub("\\]", "", rownames(x))
  rownames(x) <- gsub(" ", "_", rownames(x))
  rownames(x) <- gsub("-", "_", rownames(x))
  rownames(x) <- gsub("\\/", "_", rownames(x))
  return(x)
}

mb <- fix.names(mb)

# sort mb by most variable
mb.means <- rowMeans(mb)
mb.sd <- apply(mb, 1, function(x) sd(x))
mb.var <- mb.sd/abs(mb.means)

# get top 10 most variable
most.var <- names(head(sort(mb.var, decreasing = T),n =15))
most.common <- names(head(sort(mb.means, decreasing = T), n = 15))


#species.to.test <- c(most.var, most.common)
species.to.test <- c(most.var)

# subset mb to these and transpose for later analysis
mb.of.interest <- as.data.frame(t(mb[species.to.test,]))


### TODO: rewrite to make a list of each persons mb, mbpcs and dietpcs from the global/everyone dfs


list.mb<- NULL
list.mbpc <- NULL
list.dietpc <- NULL

for (i in seq_along(unique(map$UserName))){
  a <- as.character(unique(map$UserName)[i])
  ids <- droplevels(map$X.SampleID[map$UserName == a])
  
  sub.mb <- mb.of.interest[rownames(mb.of.interest) %in% ids,]
  list.mb[[i]] <- sub.mb

  sub.mbpc <- mbpc[rownames(mbpc) %in% ids,]
  list.mbpc[[i]] <- sub.mbpc
  
  sub.dietpc <- diet3daypc[rownames(diet3daypc) %in% ids,]
  list.dietpc[[i]] <- sub.dietpc
  
} 




# join the microbiome and diet dataframes and only keep days we have data for both
predictors <- lapply(1:34, function(i) merge(list.mbpc[[i]], list.dietpc[[i]], by = 0))


################################
# we care about the previous day ra for the model, so add these to the predictors
ra.predictors <- lapply(list.mb, function(x) as.data.frame(x))
ra.predictors <- lapply(ra.predictors, function(x) {
  x$Row.names <- rownames(x)
  return(x)
})

# we also care about these as the response variables
response <- lapply(list.mb, function(x) as.data.frame(x))
response <- lapply(response, function(x) {
  x$Row.names <- rownames(x)
  return(x)
})


# add the day variable to both predictors and response variables
predictors_day <- lapply(1:34, function(i) merge(map, predictors[[i]], by.x = "X.SampleID", by.y = "Row.names", all =F))
response_day <- lapply(1:34, function(i) merge(map, response[[i]], by.x = "X.SampleID", by.y = "Row.names", all = F))

# get the indicies to use for filtering
possible.days <- c(1:17) # vector to compare to 1:17 becuase there are 17 possible days
pred_days <- lapply(1:34, function(i) possible.days %in% predictors_day[[i]]$StudyDayNo) # this gives the list of days we values for

# make a list of the dataframes made from the vector
possible.days.df <- as.data.frame(possible.days, drop = F)

possible.days.list <- lapply(1:34, function(i) {cbind(possible.days.df, pred_days[[i]])})
possible.days.list <- lapply(possible.days.list, function(x){
  colnames(x) <- c("possible.days", "pred.days")
  return(x)
})
possible.days.list <- lapply(possible.days.list, function(x){
  x$pred_days <- ifelse(x$pred.days == TRUE, x$possible.days, NA)
  x$resp_days <- x$pred_days + c(diff(x$pred_days), NA)
  x$pred_days <- x$resp_days - 1
  keeps <- c("possible.days", "pred_days", "resp_days")
  x <- x[keeps]
  return(x)
})

# possible.days.list can now act as a framework to add the appropriate prediction and response days

# merge to get all the information in one table that we will use for testing
# either boosted regression or randomforests
dat <- lapply(1:34, function(i){merge(possible.days.list[[i]], predictors_day[[i]], by.x = "pred_days", by.y = "StudyDayNo")})
dat <- lapply(1:34, function(i){merge(dat[[i]], ra.predictors[[i]], by.x = "X.SampleID", by.y = "Row.names")})
dat <- lapply(1:34, function(i){merge(dat[[i]], response_day[[i]], by.x = "resp_days", by.y = "StudyDayNo")})

# addnames
mynames <- as.character(unique(map$UserName))

dat <- setNames(dat, mynames)

# identify people with < 5 time points
drops <- which(lapply(dat, function(x) dim(x)[1])<=5)

# drop the people with too few time points
dat <- dat[!names(dat) %in% names(drops)]

# 
# ###### Scale up for each person #########
# # read in each microbiome distance matrix and calculate pcoas for each person
# setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_1day/")
# 
# # get the file paths for the euclidean distance matrices for microbiome
# temp <- list.files(pattern = "*_tax.txt", recursive = T)
# temp <- temp[grep("euclidean", temp)]
# 
# # drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
# temp <- temp[grep("MCTs05", temp, invert = T)]
# temp <- temp[grep("MCTs06", temp, invert = T)]
# temp <- temp[grep("MCTs28", temp, invert = T)]
# temp <- temp[grep("MCTs29", temp, invert = T)]
# 
# # create a list containing each of these
# mb.dist.1day <- lapply(temp, function(x) read.delim(x, row.names = 1))
# 
# # make the PCOA from the distances
# mb.pcoa.1day <- lapply(mb.dist.1day, function(x) as.data.frame(pcoa(x)$vectors))
# 
# # limit to the top 10 axis
# mb.pcoa.1day <- lapply(mb.pcoa.1day, function (x) x[,1:10])
# mb.pcoa.1day <- lapply(mb.pcoa.1day, setNames, paste0("Mb.Axis.", 1:10))
# 
# 
# ##############################
# 
# # read in each dietary distance matrix and calculate pcoas for each person
# setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_3day/")
# 
# # get the file paths for the unweighted distance matrices for diet
# temp <- list.files(pattern = "*_food.txt", recursive = T)
# temp <- temp[grep("unweighted", temp)]
# 
# # drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
# temp <- temp[grep("MCTs05", temp, invert = T)]
# temp <- temp[grep("MCTs06", temp, invert = T)]
# temp <- temp[grep("MCTs28", temp, invert = T)]
# temp <- temp[grep("MCTs29", temp, invert = T)]
# 
# # create a list containing each of these
# diet.dist.3day <- lapply(temp, function(x) read.delim(x, row.names = 1))
# 
# # make the PCOA from the distances
# diet.pcoa.3day <- lapply(diet.dist.3day, function(x) as.data.frame(pcoa(x)$vectors))
# 
# # limit to the top 5 axis
# diet.pcoa.3day <- lapply(diet.pcoa.3day, function (x) x[,1:5])
# diet.pcoa.3day <- lapply(diet.pcoa.3day, setNames, paste0("Dt.Axis.", 1:5))
# 
# 
# # join the microbiome and diet dataframes and only keep days we have data for both
# predictors <- lapply(1:28, function(i) merge(mb.pcoa.1day[[i]], diet.pcoa.3day[[i]], by = 0))
# 
# ############################
# 
# # load the microbiome relative abundance table for each person
# # read in each dietary distance matrix and calculate pcoas for each person
# setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_1day/")
# 
# # get the file paths for the unweighted distance matrices for diet
# temp <- list.files(pattern = "*_tax.txt")
# 
# # drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
# temp <- temp[grep("MCTs05", temp, invert = T)]
# temp <- temp[grep("MCTs06", temp, invert = T)]
# temp <- temp[grep("MCTs28", temp, invert = T)]
# temp <- temp[grep("MCTs29", temp, invert = T)]
# 
# 
# # create a list containing each of these taxonomy tables
# mb.ra <- lapply(temp, function(x) read.delim(x, row.names = 1))
# 
# # fix the naming of each species
# mb.ra <- lapply(mb.ra, function(x) {
#   rownames(x) <- gsub("?.*s__", "", rownames(x))
#   rownames(x) <- gsub("?.*g__", "Uncl. Genus ", rownames(x))
#   rownames(x) <- gsub("?.*f__", "Uncl. Family ", rownames(x))
#   rownames(x) <- gsub("?.*o__", "Uncl. Order ", rownames(x))
#   rownames(x) <- gsub("?.*p__", "Uncl. Phylum ", rownames(x))
#   rownames(x) <- gsub("?.*k__", "Uncl. Kingdom ", rownames(x))
#   rownames(x) <- gsub(";NA", "", rownames(x))
#   rownames(x) <- gsub("\\[", "", rownames(x))
#   rownames(x) <- gsub("\\]", "", rownames(x))
#   rownames(x) <- gsub(" ", "_", rownames(x))
#   rownames(x) <- gsub("-", "_", rownames(x))
#   rownames(x) <- gsub("\\/", "_", rownames(x))
#   return(x)
# })
# 
# # sort each dataframe within the list by most abundant
# # sort mb.ra by most abundant
# mb.means <- lapply(mb.ra, function(x) rowMeans(x))
# 
# # get names of top 20 most abundant
# # kind of a round-about way to do this
# # but it works when pulling in each person's pre-processed taxa table
# # otherwise, have to re-code to reference the master taxa table for everyone
# most.abund <- lapply(mb.means, function(x) names(head(sort(x, decreasing = T), n = 130))) # this results in 20 common species across all 28 people
# 
# # subset mb.ra to these
# mb.ra <- lapply(1:28, function(i) mb.ra[[i]][most.abund[[i]],])
# 
# 
# #################################
# 
# # subset to just the samples for this person
# mb.ra <- lapply(1:28, function(i) mb.ra[[i]][,colnames(mb.ra[[i]]) %in% predictors[[i]]$Row.names])
# 
# 
# ################################
# # we care about the previous day ra for the model, so add these to the predictors
# ra.predictors <- lapply(mb.ra, function(x) as.data.frame(t(x)))
# ra.predictors <- lapply(ra.predictors, function(x) {
#   x$Row.names <- rownames(x)
#   return(x)
# })
# 
# # we also care about these as the response variables
# response <- lapply(mb.ra, function(x) as.data.frame(t(x)))
# response <- lapply(response, function(x) {
#   x$Row.names <- rownames(x)
#   return(x)
# })
# 
# 
# # add the day variable to both predictors and response variables
# predictors_day <- lapply(1:28, function(i) merge(map, predictors[[i]], by.x = "X.SampleID", by.y = "Row.names", all =F))
# response_day <- lapply(1:28, function(i) merge(map, response[[i]], by.x = "X.SampleID", by.y = "Row.names", all = F))
# 
# # get the indicies to use for filtering
# possible.days <- c(1:17) # vector to compare to 1:17 becuase there are 17 possible days
# pred_days <- lapply(1:28, function(i) possible.days %in% predictors_day[[i]]$StudyDayNo) # this gives the list of days we values for
# 
# # make a list of the dataframes made from the vector
# possible.days.df <- as.data.frame(possible.days, drop = F)
# 
# possible.days.list <- lapply(1:28, function(i) {cbind(possible.days.df, pred_days[[i]])})
# possible.days.list <- lapply(possible.days.list, function(x){
#   colnames(x) <- c("possible.days", "pred.days")
#   return(x)
# })
# possible.days.list <- lapply(possible.days.list, function(x){
#   x$pred_days <- ifelse(x$pred.days == TRUE, x$possible.days, NA)
#   x$resp_days <- x$pred_days + c(diff(x$pred_days), NA)
#   x$pred_days <- x$resp_days - 1
#   keeps <- c("possible.days", "pred_days", "resp_days")
#   x <- x[keeps]
#   return(x)
# })
# 
# # possible.days.list can now act as a framework to add the appropriate prediction and response days
# 
# # merge to get all the information in one table that we will use for testing
# # either boosted regression or randomforests
# dat <- lapply(1:28, function(i){merge(possible.days.list[[i]], predictors_day[[i]], by.x = "pred_days", by.y = "StudyDayNo")})
# dat <- lapply(1:28, function(i){merge(dat[[i]], ra.predictors[[i]], by.x = "X.SampleID", by.y = "Row.names")})
# dat <- lapply(1:28, function(i){merge(dat[[i]], response_day[[i]], by.x = "resp_days", by.y = "StudyDayNo")})
# 
# 

################### TEST DIET IMPACT ###########

# 
# 
# # want to look at only same speices across people
# # so find the species in everyonen - can probably do this more eligantly later if needed 
# all_species <- unlist(most.abund)
# common_species <- names(which(table(all_species)==28)) # uses common species found in everyone
# 

species <- species.to.test



"test.diet.impact" <- function(){
  set.seed(42)
  
  Master.own.diet <- NULL
  Master.others.diet <- NULL
  Master.own.microbiome <- NULL
  Master.adj.own.diet.pvalues <- NULL
  Master.diet.pc1 <- NULL
  Master.diet.pc2 <- NULL
  Master.diet.pc3 <- NULL
  Master.diet.pc1.pval <- NULL
  Master.diet.pc2.pval <- NULL
  Master.diet.pc3.pval <- NULL
  
  
  for (k in 1:length(dat)) {
    # person A index
    a <- k
    
    # person B index
    possible.b <- 1:length(dat)
    possible.b <- possible.b[possible.b != a]
    b <- sample(possible.b,1)
    
    # get the species names
    today.species <- paste0(species, ".x")
    tomorrow.species <- paste0(species, ".y")
    
    # create empty list for storing values
    own.diet.predictions <- NULL
    others.diet.predictions <- NULL
    own.microbiome.predictions <- NULL
    own.diet.pvalues <- NULL
    own.diet.pc1 <- NULL
    own.diet.pc2 <- NULL
    own.diet.pc3 <- NULL
    own.diet.pc1.pval <- NULL
    own.diet.pc2.pval <- NULL
    own.diet.pc3.pval <- NULL
    
    for (i in seq_along(tomorrow.species)) {
      
      tomorrow <- tomorrow.species[i]
      today <- today.species[i]
      
      mb.features <- c(today, "Mb.Axis.1","Mb.Axis.2", "Mb.Axis.3")
      dt.features <- c("Dt.Axis.1","Dt.Axis.2", "Dt.Axis.3")
      
      # predict person A's microbiome from microbiome alone
      # regress tomorrow (yA) from today's species and microbiome
      # predict tomorrow (yAhat)
      yA.feats <- dat[[a]][,mb.features]
      yA <- dat[[a]][,tomorrow]
      yA.model <- lm(yA ~., yA.feats)
      yAhat <- predict(yA.model)
      
      # get residuals
      yAr <- yAhat - yA
      
      # regress the residuals on diet
      modelA.feats <- dat[[a]][,dt.features]
      modelA <- lm(yAr ~., modelA.feats)
      yrAhat <- predict(modelA)
      
      # pull pcs
      own.diet.pc1[[tomorrow]] <- summary(modelA)[[4]][2,1]
      own.diet.pc2[[tomorrow]] <- summary(modelA)[[4]][3,1]
      own.diet.pc3[[tomorrow]] <- summary(modelA)[[4]][4,1]
      
      # pull pvales for pcs
      own.diet.pc1.pval[[tomorrow]] <- summary(modelA)[[4]][2,4]
      own.diet.pc2.pval[[tomorrow]] <- summary(modelA)[[4]][3,4]
      own.diet.pc3.pval[[tomorrow]] <- summary(modelA)[[4]][4,4]
      
      # subtract predicted residuals (modelA) from yhat
      # and correlate with y
      own.diet <- cor(yAhat - yrAhat, yA)
     
      
      # compare yA and yAhat to 
      own.microbiome <- cor(yAhat, yA)
      
      # delta between own diet and own mb
      own.diet.delta <- own.diet - own.microbiome
      
      ##### permutation part #####
      # predict for A from scrambled diet features
      # scramble the diet features
      scramble.cor <- NULL
      scramble.cor.delta <- NULL
      
      for (j in 1:10000) { # TODO dont forget to change back
        scramble.Afeats <- modelA.feats[sample(nrow(modelA.feats)),]
        scramble.Amodel <- lm(yAr ~., scramble.Afeats)
        yrAhat.scramble <- predict(scramble.Amodel)
        scramble.cor[[j]] <- cor(yAhat - yrAhat.scramble, yA)
        scramble.cor.delta[[j]] <- scramble.cor[[j]] - own.microbiome
      } 
      
      #own.diet.p <- sum(scramble.cor > own.diet)/10
      own.diet.p <- mean(scramble.cor.delta > own.diet.delta)
      
      
      ##
      # predict person B
      # predict person B's species from B microbiome alone
      yB.feats <- dat[[b]][,mb.features]
      yB <- dat[[b]][,tomorrow]
      yB.model <- lm(yB ~., yB.feats)
      yBhat <- predict(yB.model)
      # get residuals
      yBr <- yBhat - yB
      
      #run predict on the model with person A, but passing in the diet data for person B 
      #(you can pass in a new data matrix to predict) 
      #to get their externally-predicted residuals yr.B.hat.
      modelB.feats <- dat[[b]][,dt.features]
      yrBhat <- predict(modelA, modelB.feats)

      #subtract yr.B.hat from yhat.B and correlate the final result with y.B.
      others.diet <- cor(yBhat - yrBhat, yB)
      
     
      
      # store variables
      
      own.diet.predictions[[tomorrow]] <- own.diet
      own.diet.pvalues[[tomorrow]] <- own.diet.p
      others.diet.predictions[[tomorrow]] <- others.diet
      own.microbiome.predictions[[tomorrow]] <- own.microbiome

    }
    Master.own.diet[[k]] <- own.diet.predictions
    Master.others.diet[[k]] <- others.diet.predictions
    Master.own.microbiome[[k]] <- own.microbiome.predictions
    Master.diet.pc1[[k]] <- own.diet.pc1
    Master.diet.pc2[[k]] <- own.diet.pc2
    Master.diet.pc3[[k]] <- own.diet.pc3
    Master.diet.pc1.pval[[k]] <- own.diet.pc1.pval
    Master.diet.pc2.pval[[k]] <- own.diet.pc2.pval
    Master.diet.pc3.pval[[k]] <- own.diet.pc3.pval
    
    
    # TODO: will only adjust when rerun with perm = 10000
    # for now just running with 99 perms for debugging
    
    Master.adj.own.diet.pvalues[[k]] <- p.adjust(own.diet.pvalues, "fdr")  

  }
  
  return(list(Master.own.diet = Master.own.diet, 
              Master.others.diet = Master.others.diet, 
              Master.own.microbiome = Master.own.microbiome, 
              Master.adj.own.diet.pvalues = Master.adj.own.diet.pvalues,
              Master.diet.pc1 = Master.diet.pc1,
              Master.diet.pc2 = Master.diet.pc2,
              Master.diet.pc3 = Master.diet.pc3,
              Master.diet.pc1.pval = Master.diet.pc1.pval,
              Master.diet.pc2.pval = Master.diet.pc2.pval,
              Master.diet.pc3.pval = Master.diet.pc3.pval))
  
}


myresults <- test.diet.impact()
save(myresults, file ='output/permutation.results.var.Rdata')

#load(file = "../../../../output/permutation.results.Rdata")

# quick look at the means
mean.own.diet <- unlist(lapply(myresults$Master.own.diet, function(x) mean(x)))
mean.others.diet <- unlist(lapply(myresults$Master.others.diet, function(x) mean(x)))
mean.own.microbiome <- unlist(lapply(myresults$Master.own.microbiome, function(x) mean(x)))


# format myresults as dataframes for plotting
own.diet <- data.frame(matrix(unlist(myresults$Master.own.diet), nrow = length(species), byrow = F))
rownames(own.diet) <- species


mynames <- gsub("MCTs", "", names(dat))
colnames(own.diet) <- mynames

others.diet <- data.frame(matrix(unlist(myresults$Master.others.diet), nrow = length(species), byrow = F))
rownames(others.diet) <- species
colnames(others.diet) <- mynames

own.microbiome <- data.frame(matrix(unlist(myresults$Master.own.microbiome), nrow = length(species), byrow = F))
rownames(own.microbiome) <- species
colnames(own.microbiome) <- mynames

own.diet.pvalues <- data.frame(matrix(unlist(myresults$Master.adj.own.diet.pvalues), nrow = length(species), byrow = F))
rownames(own.diet.pvalues) <- species
colnames(own.diet.pvalues) <- mynames


own.diet.pc1 <- data.frame(matrix(unlist(myresults$Master.diet.pc1), nrow = length(species), byrow = F))
rownames(own.diet.pc1) <- species
colnames(own.diet.pc1) <- mynames

own.diet.pc2 <- data.frame(matrix(unlist(myresults$Master.diet.pc2), nrow = length(species), byrow = F))
rownames(own.diet.pc2) <- species
colnames(own.diet.pc2) <- mynames

own.diet.pc3 <- data.frame(matrix(unlist(myresults$Master.diet.pc3), nrow = length(species), byrow = F))
rownames(own.diet.pc3) <- species
colnames(own.diet.pc3) <- mynames


own.diet.pc1.pvals <- data.frame(matrix(unlist(myresults$Master.diet.pc1.pval), nrow = length(species), byrow = F))
rownames(own.diet.pc1.pvals) <- species
colnames(own.diet.pc1.pvals) <- mynames


own.diet.pc2.pvals <- data.frame(matrix(unlist(myresults$Master.diet.pc2.pval), nrow = length(species), byrow = F))
rownames(own.diet.pc2.pvals) <- species
colnames(own.diet.pc2.pvals) <- mynames


own.diet.pc3.pvals <- data.frame(matrix(unlist(myresults$Master.diet.pc3.pval), nrow = length(species), byrow = F))
rownames(own.diet.pc3.pvals) <- species
colnames(own.diet.pc3.pvals) <- mynames






# format for ggplot to make plots
require(ggplot2)
require(reshape2)

own.diet.melt <- rownames_to_column(own.diet, var = "Taxa")
own.diet.melt <- melt(own.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.corr")

others.diet.melt <- rownames_to_column(others.diet, var = "Taxa")
others.diet.melt <- melt(others.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "others.diet.corr")

own.microbiome.melt <-rownames_to_column(own.microbiome, var = "Taxa")
own.microbiome.melt <- melt(own.microbiome.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.mb.corr")

own.diet.pvalues <- rownames_to_column(own.diet.pvalues, var = "Taxa")
own.diet.pvalues.melt <- melt(own.diet.pvalues, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.p")

own.diet.pc1 <- rownames_to_column(own.diet.pc1, var = "Taxa")
own.diet.pc1.melt <- melt(own.diet.pc1, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.pc1")

own.diet.pc2 <- rownames_to_column(own.diet.pc2, var = "Taxa")
own.diet.pc2.melt <- melt(own.diet.pc2, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.pc2")

own.diet.pc3 <- rownames_to_column(own.diet.pc3, var = "Taxa")
own.diet.pc3.melt <- melt(own.diet.pc3, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.pc3")

own.diet.pc1.pvals <- rownames_to_column(own.diet.pc1.pvals, var = "Taxa")
own.diet.pc1.pvals.melt <- melt(own.diet.pc1.pvals, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.pc1.pvals")

own.diet.pc2.pvals <- rownames_to_column(own.diet.pc2.pvals, var = "Taxa")
own.diet.pc2.pvals.melt <- melt(own.diet.pc2.pvals, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.pc2.pvals")

own.diet.pc3.pvals <- rownames_to_column(own.diet.pc3.pvals, var = "Taxa")
own.diet.pc3.pvals.melt <- melt(own.diet.pc3.pvals, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.pc3.pvals")


# make master df

master.df <- merge(own.diet.melt,others.diet.melt)
master.df <- merge(own.microbiome.melt, master.df)

master.df.melt <- melt(master.df, id.vars = c("Taxa", "UserName"))

# get summary stats for plotting error bars
summary.master <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "mean")
colnames(summary.master) <- c("variable", "mean")
summary.sd <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "sd")
colnames(summary.sd) <- c("variable", "sd")

summary.master <- merge(summary.master, summary.sd)
summary.master$se <- summary.master$sd/sqrt(28) # se is sd/sqrt(n)

# make the plot
mybarplot <- ggplot(data = master.df.melt, aes(x = variable, y = value)) +
  geom_bar(stat = "summary", fun.y = mean, aes(fill = as.factor(variable))) +
  geom_errorbar(data = summary.master, aes(x = variable, ymin = mean - se, ymax = mean + se), width = 0.4, inherit.aes = F) +
  theme_classic()

# pring the plot
mybarplot



# make paired difference plots by person

require(gridExtra)
require(directlabels)

mb.dt <- master.df.melt[master.df.melt$variable %in% c("own.mb.corr", "own.diet.corr"),]

mb.dt.wide <- merge(master.df, own.diet.pvalues.melt)
mb.dt.wide <- mb.dt.wide[,!colnames(mb.dt.wide) %in% c("others.diet.corr")]
mb.dt.wide$sig<- ifelse(mb.dt.wide$own.diet.p <= 0.25, "sig", "ns")
mb.dt.wide$groups <- ifelse(mb.dt.wide$sig == "sig", mb.dt.wide$Taxa, NA)

mb.dt.for.merge <- mb.dt.wide[colnames(mb.dt.wide) %in% c("Taxa", "UserName", "groups")]

mb.dt <- merge(mb.dt, mb.dt.for.merge)


# make multiple plots in a loop

single_plots <- as.list(NULL)

plot <- mb.dt[order(mb.dt$UserName, decreasing = F),]
plot <- plot[order(plot$groups, na.last = F),]

group.names <- unique(mb.dt$groups)
group.colors <- c("#858c69", "#ff40f2", "#ff0000", "#008c4b", "#00138c","#8c696e", "#ffbfbf", "#8c7723", "#468c75", "#8091ff",
                  "#ff80c4", "#8c3123", "#fff2bf", "#40fff2", "#69698c", "#ff0044", "#ff9180", "#e5ff80", "#bffbff", "#5940ff",
                  "#858c69", "#40d9ff", "#c480ff", "#ff8c40", "#4b8c00", "#23698c", "#69238c", "#8c4b00", "#bfffbf", "#004b8c",
                  "#eabfff")# "#ffc480", "#40ff59", "#80c4ff", "#ffd940" )
names(group.colors) <- group.names


# view colors
cols <- function(a) image(1:31, 1, as.matrix(1:31), col=a, axes=FALSE , xlab="", ylab="")
cols(group.colors)

for (i in 1:length(unique(plot$UserName))) {
  j <- as.character(unique(plot$UserName)[[i]])
  plot_un <- plot[plot$UserName == j,]
  
  plot_un <- plot_un[order(plot_un$groups, na.last = T),]
  
  single_plots[[i]] <- ggplot(data = plot_un, aes(x = variable, y = value)) + 
    geom_line(aes(group = Taxa), color = "grey", alpha = 0.5) + 
    geom_line(aes(group = Taxa, color = factor(groups)), size = 1) +
    scale_color_manual(values = group.colors) +
    #geom_dl(aes(label = groups), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.5, hjust = 1)) +
    geom_dl(aes(label = groups), method = "last.polygons") +
    geom_point() +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(size = 7),
          legend.title = element_blank(),
          axis.title = element_blank()) +
    annotate("text", x = 1.5, y = 0, label = j)
  
}


grid.arrange(grobs = single_plots[1:28], nrow=7)




#### make plots of pcs and pvals by species
require(ggbeeswarm)
modelpcs <- merge(own.diet.pc1.melt, own.diet.pc1.pvals.melt)
modelpcs <- merge(modelpcs, own.diet.pc2.melt)
modelpcs <- merge(modelpcs, own.diet.pc2.pvals.melt)
modelpcs <- merge(modelpcs, own.diet.pc3.melt)
modelpcs <- merge(modelpcs, own.diet.pc3.pvals.melt)

modelpcs$pc1sig <- ifelse(modelpcs$own.diet.pc1.pvals < 0.05, "sig", "ns")
modelpcs$pc2sig <- ifelse(modelpcs$own.diet.pc2.pvals < 0.05, "sig", "ns")
modelpcs$pc3sig <- ifelse(modelpcs$own.diet.pc3.pvals < 0.05, "sig", "ns")

# reorder taxa for plotting
modelpcs$Taxa <- factor(modelpcs$Taxa, levels = species.to.test)

pcplot <- ggplot(data = modelpcs, aes(x = Taxa, y = own.diet.pc1)) +
  #geom_violin(trim = F) +
  geom_quasirandom(data = modelpcs, aes(color = pc1sig), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("grey", "blue")) +
  geom_label(data = subset(modelpcs, pc1sig == "sig"), aes(label = UserName), size = 2, nudge_y = 7, nudge_x = 0.3) +
  ylim(-150,150) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "grey"))

pcplot

# 
# pcplot2 <- ggplot(data = modelpcs, aes(x = Taxa, y = own.diet.pc2)) +
#   geom_violin(trim = F) +
#   geom_quasirandom(data = modelpcs, aes(color = pc2sig), size = 2, alpha = 0.8) +
#   geom_label(data = subset(modelpcs, pc2sig == "sig"), aes(label = UserName), size = 2, nudge_y = 7, nudge_x = 0.3) +
#   ylim(-150,150) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.major = element_line(colour = "grey"))
# 
# pcplot2
# 
# pcplot3 <- ggplot(data = modelpcs, aes(x = Taxa, y = own.diet.pc3)) +
#   geom_violin(trim = F) +
#   geom_quasirandom(data = modelpcs, aes(color = pc3sig), size = 2, alpha = 0.8) +
#   geom_label(data = subset(modelpcs, pc3sig == "sig"), aes(label = UserName), size = 2, nudge_y = 7, nudge_x = 0.3) +
#   ylim(-150,150) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.major = element_line(colour = "grey"))
# 
# pcplot3

## TODO: re do with one golbal PCOA for diet and mB (CHECK!)
# coef. for significance (plot coef vs species [x=species, y = coef and color by significance], 
# cluster > 0 then Diet microbe relationship, cluster = 0, no relationship, two clusters = personalized)
# cocurrent look at MC with new axis (CHECK, rerun with propper permutations)
# try with 2 and try with 3 for both diet and mb to avoid overfitting
# do the soylent people behave differently?

# cant do formal personalized microbiome in a formal way. Nothing conclusive because we don't have people on the same diets
# However we do have these two soylent subjects and we can show something for them
# show bar plots for the two soylent people


