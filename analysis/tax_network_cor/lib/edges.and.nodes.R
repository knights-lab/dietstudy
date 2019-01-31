### Updating to scale up from 1 person's network to the average network for the study!

require(reshape2)

setwd("/Users/abby/Documents/Projects/dietstudy/analysis/tax_network_cor/sparcc_results/")

### Import correlation matrices ###
tmp1 <- list.files(pattern = "*_out.txt")
sparcc <- lapply(tmp1, function(x) read.delim(x, row.names = 1))

### Import pvalue matrices ###
tmp2 <- list.files(pattern = "*_two_sided.txt")
pvals <- lapply(tmp2, function(x) read.delim(x, row.names = 1))

## ONE 
# sparcc <- as.matrix(sparcc)
# sparcc[lower.tri(sparcc)] <- NA
# sparcc <- melt(sparcc, na.rm = TRUE)

# pvals <- as.matrix(pvals)
# pvals[lower.tri(pvals)] <- NA
# pvals <- melt(pvals)
# pvals <- na.omit(pvals)

# sparcc <- cbind(pvals, sparcc[3])
# colnames(sparcc) <- c("Source", "Target", "Pval", "Sparcc" )

## MANY
sparcc <- lapply(sparcc, as.matrix)                                  # convert to matrix
sparcc <- lapply(sparcc, function(xx) {xx[lower.tri(xx)] <- NA; xx}) # set lower tri to NA
sparcc <- lapply(sparcc, function(xxi) {melt(xxi, na.rm = TRUE)})    # make long form and drop NA

pvals <- lapply(pvals, as.matrix)                                    # convert to matrix
pvals <- lapply(pvals, function(xx) {xx[lower.tri(xx)] <- NA; xx})   # set lower tri to NA
pvals <- lapply(pvals, function(xxi) {melt(xxi, na.rm = TRUE)})      # make long form and drop NA
pvals <- lapply(pvals, function(xi) {xi[,3]})                        # limit to just pvals


sparcc <- mapply(cbind, sparcc, pvals, SIMPLIFY = FALSE)             # add pvals to correlations
sparcc <- lapply(sparcc, function(x) {colnames(x) <- c("Source", "Target", "Sparcc", "Pval"); x})


# Process for averaging
# 1. Set non-sig correlations to 0
# If pvalue is significant, new = sparcc as is, else new = 0
allnew <- lapply(sparcc, function(x) {transform (x, new = ifelse(x[,"Pval"] <= 0.25, x[,"Sparcc"], 0))})

# 2. Average correlations and pvals
allnew <- lapply(allnew, function(xi) {xi[,4:5]}) 
allnew <- abind(allnew, along = 3)
allnew <- apply(allnew, c(1,2), function(x) mean(x, na.rm = TRUE))
allnew <- as.data.frame(allnew)
allnew$Source <- sparcc[[1]]$Source # we haven't removed anything, so names are the same
allnew$Target <- sparcc[[2]]$Target # we haven't removed anything, so names are the same

# 3. limit to ave pvals <=0.25, all correlations are tiny!
allnew <- allnew[allnew$Pval <= 0.17,]

### Make the edge table ###
edges <- allnew
nodes <- c(as.character(edges$Source), as.character(edges$Target))
nodes <- unique(nodes)
nodes <- as.data.frame(nodes)
nodes$taxonomy <- nodes$nodes
nodes$taxonomy <- gsub("\\.", ";", nodes$taxonomy)
nodes <- nodes %>% separate(col = taxonomy, into = c("L1", "L2", "L3", "L4", "L5", "L6"), sep = ";")


# fix naming for edges and nodes
edges$Source <- gsub(".*s__", "", edges$Source)
edges$Target <- gsub(".*s__", "", edges$Target)
nodes$nodes <- gsub(".*s__", "", nodes$nodes)


write.table(edges, file = "../ave.edge.csv", sep = ",", quote = F, row.names = F)
write.table(nodes, file = "../ave.nodes.csv", sep = ",", quote = F, row.names = F)


