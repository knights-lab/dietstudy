require(dplyr)
require(ggplot2)

setwd("/Users/abby/Documents/Projects/dietstudy/")

map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
map_username <- read.table("data/maps/UserName_map.txt", sep = "\t", header = T, comment = "")

### create a significant procrustes variable ###
procrustes.pvals <- read.table("analysis/procrustes/results_username/totals.txt", comment = "#", header = F, sep = "\t")
colnames(procrustes.pvals) <- c("FP1","FP2","Num included dimensions", "Monte Carlo p-value", "Count better", "M^2")
procrustes.pvals$UserName <- gsub("\\_.*", "", procrustes.pvals$FP1)
procrustes.pvals <- procrustes.pvals %>% select(UserName, `Monte Carlo p-value`, `Count better`, `M^2`)
procrustes.pvals$Procrustes <- ifelse(procrustes.pvals$`Monte Carlo p-value` <= 0.06, "sig", "ns")

map_username <- full_join(map_username, procrustes.pvals)

### create variable to indicate confidence in dietary reporting (high/low) ###
  
  # load nutr data
  nutr <- read.table("data/processed_nutr/nutr_totals.txt", sep = "\t", header = T, comment = "")

  # get average kcal intake per person
  mean_kcal <- nutr %>% group_by(UserName) %>% summarize(mean_kcal = mean(KCAL))
  
  # join with map_username
  map_username <- full_join(map_username, mean_kcal, by = "UserName")
  
  # calculate REE and TEE
  map_username <- map_username %>% mutate(REE = ifelse(Gender=="M", (10*Weight+6.25*Height-5*Age+5), 
                                                       ifelse(Gender=="F", (10*Weight+6.25*Height-5*Age-161),
                                                              NA)), 
                                          TEE = REE*Activity.Factor) # calculate TEE with activity factor

  # Energy intake difference actual - estimated 
  map_username$kcal_diff <-  map_username$mean_kcal - map_username$TEE 
  

### Write new mapping file ###  
write.table(map_username, "data/maps/UserName_map.txt", sep = "\t", col.names = T, row.names = F, quote = F)




