strains = read.delim("mwa.str",row=1); # mwa is embamulate output from burst "otu"
map = read.delim('SampleID_map.txt',row=1,as.is=T);

# Manicure the samplenames, grab latest (mct study only!)
colnames(strains) = gsub(".S[0-9]+.R1.001","",colnames(strains));      # Clean old plate IDs
strains = strains[,order(colnames(strains))];              # Sort nicely by sample ID
strains = strains[,-(grep("L0",colnames(strains))-1)];     # Keep new runs only
colnames(strains) = gsub(".S[0-9]+.L001.R1.001","",colnames(strains)); # Clean new plate IDs

# Filter samps and strains
map = map[map$StudyDayNo >= 10 | map$StudyDayNo <= 7,]       # chew off the supplement transition days
map = map[map$Study.Status != "Dropped",]
#strains = strains[,colSums(strains) > 20000]
matched = intersect(colnames(strains),rownames(map))         # get samples in common (table & map)
map = map[matched,]; strains = as.matrix(strains[,matched])  # Do the match for both
strains = strains[rowSums(strains) >= sum(strains)/1000000 &
                    rowSums(strains>0) > 5,,drop=F] # Keep strains in >= a millionth of the data
strains.c = strains

# Transform abundances
strains = sweep(strains,2,colSums(strains),'/');
strains = sweep(sqrt(strains),2,colSums(sqrt(strains)),'/');

# OR do CLR!
strains = t(strains); eps = 0.5
strains = strains*(1 - rowSums(strains==0)*eps/rowSums(strains))
strains[strains==0]=eps
strains = sweep(strains,1,rowSums(strains),'/');
ls = log(strains)
strains = t(ls - rowMeans(ls))
strains = strains[,!is.nan(colSums(strains))]

# Aggregate map and taxa table by username 
map.ag = aggregate(map, by=list(map$UserName, map$StudyDayNo >= 10), FUN=min) # FUN=function(xx) min(as.character(xx)))  #function(xx) xx[which(!duplicated(xx))])
taxa.ag = aggregate(t(strains), by=list(map$UserName, map$StudyDayNo >= 10), FUN=mean)
#rownames(map.ag) = map.ag$Group.1; map.ag = map.ag[-1]
#rownames(taxa.ag) = taxa.ag$Group.1; taxa.ag = t(data.matrix(taxa.ag[,-1]))
map = map.ag
strains = t(data.matrix(taxa.ag[,3:ncol(taxa.ag)])) # col names are meaningless

# Variance feature selection
##CV = apply(strains,1,function(x) sd(x)/mean(x))
#V =  apply(strains,1,var)
#select = V >= summary(V)["Median"] #"3rd Qu."]
#strains = strains[select,]
#strains = sweep(strains,2,colSums(strains),'/');
