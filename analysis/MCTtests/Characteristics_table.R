library(plyr)
library(dplyr)
library(broom)


setwd("/Users/abby/Documents/Projects/MCTs/")

map <- read.table("Publication/Data/map.csv", sep = ",", header = TRUE)
blood_pre <- read.table("Data/Blood_draw_02-06-17_Knights.txt", sep = "\t", header = TRUE, comment = "")
blood_post <- read.table("Data/Blood_draw_02-17-17_Knights.txt", sep = "\t", header = TRUE, comment = "")

all_blood <- merge(blood_pre, blood_post)
map_all <- merge(map, all_blood)

smry <- ddply(map_all, .(Gender, Supplement), summarise, 
              Age_mean=mean(Age), 
              Age_sd = sd(Age), 
              Gender_count=length(Gender),
              Weight_mean = mean(Weight),
              Weight_sd = sd(Weight),
              Height_mean = mean(Height),
              Height_sd = sd(Height),
              Waist_cir_mean = mean(Waist.Circumference, na.rm = TRUE),
              Waist_cir_sd = sd(Waist.Circumference, na.rm = TRUE),
              Base_Chole_mean = mean(Cholesterol.Baseline, na.rm = TRUE),
              Base_Chole_sd = sd(Cholesterol.Baseline, na.rm = TRUE),
              Base_Trig_mean = mean(Trigs.Baseline, na.rm = TRUE),
              Base_Trig_sd = sd(Trigs.Baseline, na.rm = TRUE),
              Base_HDL_mean = mean(HDL.Baseline, na.rm = TRUE),
              Base_HDL_sd = sd(HDL.Baseline, na.rm = TRUE),
              Base_LDL_mean = mean(LDL.Baseline, na.rm = TRUE),
              Base_LDL_sd = sd(LDL.Baseline, na.rm = TRUE),
              Base_Glu_mean = mean(Glu.Baseline, na.rm = TRUE),
              Base_Glu_sd = sd(Glu.Baseline, na.rm = TRUE),
              Base_Ins_mean = mean(Ins.Baseline, na.rm = TRUE),
              Base_Ins_sd = sd(Ins.Baseline, na.rm = TRUE),
              Base_HOMA_IR_mean = mean(Homa.IR.Baseline, na.rm = TRUE),
              Base_HOMA_IR_sd = sd(Homa.IR.Baseline, na.rm = TRUE))

smry_supp <- ddply(map_all, .(Supplement), summarise, 
              Age_mean=mean(Age), 
              Age_sd = sd(Age), 
              Gender_count=length(Gender),
              Weight_mean = mean(Weight),
              Weight_sd = sd(Weight),
              Height_mean = mean(Height),
              Height_sd = sd(Height),
              Waist_cir_mean = mean(Waist.Circumference, na.rm = TRUE),
              Waist_cir_sd = sd(Waist.Circumference, na.rm = TRUE),
              Base_Chole_mean = mean(Cholesterol.Baseline, na.rm = TRUE),
              Base_Chole_sd = sd(Cholesterol.Baseline, na.rm = TRUE),
              Base_Trig_mean = mean(Trigs.Baseline, na.rm = TRUE),
              Base_Trig_sd = sd(Trigs.Baseline, na.rm = TRUE),
              Base_HDL_mean = mean(HDL.Baseline, na.rm = TRUE),
              Base_HDL_sd = sd(HDL.Baseline, na.rm = TRUE),
              Base_LDL_mean = mean(LDL.Baseline, na.rm = TRUE),
              Base_LDL_sd = sd(LDL.Baseline, na.rm = TRUE),
              Base_Glu_mean = mean(Glu.Baseline, na.rm = TRUE),
              Base_Glu_sd = sd(Glu.Baseline, na.rm = TRUE),
              Base_Ins_mean = mean(Ins.Baseline, na.rm = TRUE),
              Base_Ins_sd = sd(Ins.Baseline, na.rm = TRUE),
              Base_HOMA_IR_mean = mean(Homa.IR.Baseline, na.rm = TRUE),
              Base_HOMA_IR_sd = sd(Homa.IR.Baseline, na.rm = TRUE))


map_no_drop <- map_all %>% filter(Supplement == "EVOO" | Supplement == "MCT")

# Check to see if the extra level is gone
levels(map_no_drop$Supplement)

# It's not, so need to reset the factor
map_no_drop <- droplevels(map_no_drop)

levels(map_no_drop$Supplement)

# view the categories and see if there are any differences between groups at baseline
plot(Age ~ Supplement, data = map_no_drop)
t.test(Age ~ Supplement, data = map_no_drop) # p-val = 0.34

# view Gender
plot(Gender ~ Supplement, data = map_no_drop)
table(map_no_drop$Gender)

# what is the easiest way to test frequency? Chi-sq?
# TODO, figure this out??!

# view weight
plot(Weight ~ Supplement, data = map_no_drop)
t.test(Weight ~ Supplement, data = map_no_drop) # p-val = 0.79

# view height
plot(Height ~ Supplement, data = map_no_drop)
t.test(Height ~ Supplement, data = map_no_drop) # p-val = 0.87

# view waist circumference
plot(Waist.Circumference ~ Supplement, data = map_no_drop)
t.test(Waist.Circumference ~ Supplement, data = map_no_drop) # p-val = 0.50

# view baseline cholesterol
plot(Cholesterol.Baseline ~ Supplement, data = map_no_drop)
t.test(Cholesterol.Baseline ~ Supplement, data = map_no_drop) # p-val = 0.34

# view baseline Trigs
plot(Trigs.Baseline ~ Supplement, data = map_no_drop)
t.test(Trigs.Baseline ~ Supplement, data = map_no_drop) # p-val = 0.79

# view hdl baseline
plot(HDL.Baseline ~ Supplement, data = map_no_drop)
t.test(HDL.Baseline ~ Supplement, data =  map_no_drop) # p-val = 0.36

# view ldl baseline
plot(LDL.Baseline ~ Supplement, data = map_no_drop)
t.test(LDL.Baseline ~ Supplement, data = map_no_drop) # p-val = 0.49

# view glucose
plot(Glu.Baseline ~ Supplement, data = map_no_drop)
t.test(Glu.Baseline ~ Supplement, data = map_no_drop) # p-val = 0.31

# view insulin 
plot(Ins.Baseline ~ Supplement, data = map_no_drop) # p-val = 0.52
t.test(Ins.Baseline ~ Supplement, data = map_no_drop)

# At baseline, there are no differences in characteristics between supplementation groups


