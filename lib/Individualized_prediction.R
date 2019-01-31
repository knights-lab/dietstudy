# NOTE: This was working code. This has now been replaced by the test_personal_diet_plus_permutation_script


# is prediction of specific species RA from microbiome day prior improved by adding personal diet information?

# wrangling needs:
# relative abundances of the top ~25 most abundant speices
# principal coordinates from today's microbiome, first 5 PCAs
# principal coordinates from today's (3 day aggregate) diet, first 5 PCAs
# today's microbe relative abundance
# predict tomorrow's microbiome relative abundance

# example set up:
# create a dataframe containing
# B vulgatus today MB.PC1.day	MB.PC2.day	MB.PC3.day	Diet.Mean.PC1.last.3.days	Diet.Mean.PC2.last.3.days	Diet.Mean.PC3.last.3.days	B. vulgatus tomorrow
# Each row is a microbiome sample
# Do some fancy wrangling to make sure we don't include days with mismatches when there are missing fecal samples


setwd(dir = "/Users/abby/Documents/Projects/dietstudy/")
require(ape)
require(tidyverse)

# load map
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID", "StudyDayNo","Sample.Week.Day")]

###### Scale up for each person #########
# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_1day/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.1day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.1day <- lapply(mb.dist.1day, function(x) as.data.frame(pcoa(x)$vectors))

# limit to the top 10 axis
mb.pcoa.1day <- lapply(mb.pcoa.1day, function (x) x[,1:10])
mb.pcoa.1day <- lapply(mb.pcoa.1day, setNames, paste0("Mb.Axis.", 1:10))


##############################

# read in each dietary distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_3day/")

# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.3day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.3day <- lapply(diet.dist.3day, function(x) as.data.frame(pcoa(x)$vectors))

# limit to the top 5 axis
diet.pcoa.3day <- lapply(diet.pcoa.3day, function (x) x[,1:5])
diet.pcoa.3day <- lapply(diet.pcoa.3day, setNames, paste0("Dt.Axis.", 1:5))


# join the microbiome and diet dataframes and only keep days we have data for both
predictors <- lapply(1:29, function(i) merge(mb.pcoa.1day[[i]], diet.pcoa.3day[[i]], by = 0))

############################

# load the microbiome relative abundance table for each person
# read in each dietary distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy/figures/figure 4/procrustes/data_username_1day/")

# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_tax.txt")

# drop MCTs05, MCTs06, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]


# create a list containing each of these taxonomy tables
mb.ra <- lapply(temp, function(x) read.delim(x, row.names = 1))

# fix the naming of each species
mb.ra <- lapply(mb.ra, function(x) {
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
})

# sort each dataframe within the list by most abundant 
# sort mb.ra by most abundant
mb.means <- lapply(mb.ra, function(x) rowMeans(x))

# get names of top 20 most abundant
most.abund <- lapply(mb.means, function(x) names(head(sort(x, decreasing = T), n = 20)))

# subset mb.ra to these
mb.ra <- lapply(1:29, function(i) mb.ra[[i]][most.abund[[i]],])


#################################

# subset to just the samples for this person
mb.ra <- lapply(1:29, function(i) mb.ra[[i]][,colnames(mb.ra[[i]]) %in% predictors[[i]]$Row.names])


################################
# we care about the previous day ra for the model, so add these to the predictors
ra.predictors <- lapply(mb.ra, function(x) as.data.frame(t(x)))
ra.predictors <- lapply(ra.predictors, function(x) {
  x$Row.names <- rownames(x)
  return(x)
})

# we also care about these as the response variables
response <- lapply(mb.ra, function(x) as.data.frame(t(x)))
response <- lapply(response, function(x) {
  x$Row.names <- rownames(x)
  return(x)
})


# add the day variable to both predictors and response variables
predictors_day <- lapply(1:29, function(i) merge(map, predictors[[i]], by.x = "X.SampleID", by.y = "Row.names", all =F))
response_day <- lapply(1:29, function(i) merge(map, response[[i]], by.x = "X.SampleID", by.y = "Row.names", all = F))

# get the indicies to use for filtering
possible.days <- c(1:17) # vector to compare to 1:17 becuase there are 17 possible days
pred_days <- lapply(1:29, function(i) possible.days %in% predictors_day[[i]]$StudyDayNo) # this gives the list of days we values for

# make a list of the dataframes made from the vector
possible.days.df <- as.data.frame(possible.days, drop = F)

possible.days.list <- lapply(1:29, function(i) {cbind(possible.days.df, pred_days[[i]])})
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
dat <- lapply(1:29, function(i){merge(possible.days.list[[i]], predictors_day[[i]], by.x = "pred_days", by.y = "StudyDayNo")})
dat <- lapply(1:29, function(i){merge(dat[[i]], ra.predictors[[i]], by.x = "X.SampleID", by.y = "Row.names")})
dat <- lapply(1:29, function(i){merge(dat[[i]], response_day[[i]], by.x = "resp_days", by.y = "StudyDayNo")})


# so the relative abundances that end in .x are predictors, and the .y are the response days

## TODO: scale up to run on the total dat list


################ PREDICTIONS ####################


"run.prediction" <- function(){
  
  set.seed(42)
  
  master.diet.corrs <- NULL
  master.mb.corrs <- NULL
  
  for (k in 1:29) {
    
    today.species <- paste0(most.abund[[k]], ".x")
    tomorrow.species <- paste0(most.abund[[k]], ".y")
    
    diet.corrs <- NULL
    mb.corrs <- NULL
    
    
    for (i in 1:length(tomorrow.species)) {
      
      tomorrow <- tomorrow.species[i]
      today <- today.species[i]
      
      # predict each taxon with microbiome features alone
      mb.features <- c(today, "Mb.Axis.1","Mb.Axis.2")
      
      mymodel.mb <- lm(dat[[k]][,tomorrow] ~ .,data = dat[[k]][mb.features])
      predicted.mb <- predict(mymodel.mb) 
      
      # predict each taxon with first 5 microbiome and 5 diet features
      diet.features <- c(today,"Mb.Axis.1","Mb.Axis.2","Dt.Axis.1","Dt.Axis.2")
      
      mymodel.diet <- lm(dat[[k]][,tomorrow] ~ ., data = dat[[k]][diet.features])
      predicted.diet <- predict(mymodel.diet)
      
      # compare predicted to actual for mb alone
      mycorr.mb <- cor(predicted.mb, dat[[k]][,tomorrow])
      
      # compare predicted to acutal for mb + diet
      mycorr.diet <- cor(predicted.diet, dat[[k]][,tomorrow])
      
      # output the predicted v. actual correlations
      
      mb.corrs[[tomorrow]] <- mycorr.mb
      diet.corrs[[tomorrow]] <- mycorr.diet
      }
    master.mb.corrs[[k]]=mb.corrs
    master.diet.corrs[[k]]=diet.corrs
  }
  return(list(master.diet.corrs = master.diet.corrs, master.mb.corrs = master.mb.corrs))
}


myresults <- run.prediction()

diet.means <- unlist(lapply(myresults$master.diet.corrs, function(x) mean(x)))
mb.means <- unlist(lapply(myresults$master.mb.corrs, function(x) mean(x)))

boxplot(mb.means, diet.means)


