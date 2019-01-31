setwd("~/Documents/Projects/dietstudy/")

# plot predictions results
load(file = "output/permutation.results.Rdata")
load(file = "data/species.to.test.Rdata")
load(file = "data/test_personal_diet_dat.Rdata")


species <- species.to.test

# format myresults as dataframes for plotting
# Note: should really be a function - but this will do.

mynames <- gsub("MCTs", "", names(dat))


own.diet <- data.frame(matrix(unlist(myresults$Master.own.diet), nrow = length(species), byrow = F))
rownames(own.diet) <- species
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

#fdr adjust within person
own.diet.pc1.pvals <- apply(own.diet.pc1.pvals, 2, function(x) p.adjust(x, method = "fdr")) # column-wise within a person
own.diet.pc1.pvals <- as.data.frame(own.diet.pc1.pvals)

own.diet.pc2.pvals <- data.frame(matrix(unlist(myresults$Master.diet.pc2.pval), nrow = length(species), byrow = F))
rownames(own.diet.pc2.pvals) <- species
colnames(own.diet.pc2.pvals) <- mynames

#fdr adjust within person
own.diet.pc2.pvals <- apply(own.diet.pc2.pvals, 2, function(x) p.adjust(x, method = "fdr")) # column-wise within a person
own.diet.pc2.pvals <- as.data.frame(own.diet.pc2.pvals)

own.diet.pc3.pvals <- data.frame(matrix(unlist(myresults$Master.diet.pc3.pval), nrow = length(species), byrow = F))
rownames(own.diet.pc3.pvals) <- species
colnames(own.diet.pc3.pvals) <- mynames


#fdr adjust within person
own.diet.pc3.pvals <- apply(own.diet.pc3.pvals, 2, function(x) p.adjust(x, method = "fdr")) # column-wise within a person
own.diet.pc3.pvals <- as.data.frame(own.diet.pc3.pvals)


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

require(agricolae)
#Tukey
hsd=HSD.test(aov(value~variable,data=master.df.melt), "variable", group=T)

# get summary stats for plotting error bars
summary.master <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "mean")
colnames(summary.master) <- c("variable", "mean")
summary.sd <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "sd")
colnames(summary.sd) <- c("variable", "sd")

summary.master <- merge(summary.master, summary.sd)
summary.master$se <- summary.master$sd/sqrt(28) # se is sd/sqrt(n)
summary.master <- merge(summary.master, as.data.frame(hsd$groups), by.x = "mean", by.y = "value")

summary.master

# make the plot
mybarplot <- ggplot(data = master.df.melt, aes(x = factor(variable, levels = c("own.diet.corr", "own.mb.corr", "others.diet.corr")), y = value)) +
  geom_bar(stat = "summary", fun.y = mean, aes(fill = as.factor(variable)), show.legend = F) +
  scale_x_discrete(labels = c( "Diet and\nmicrobiome","Microbiome\nonly", "Microbiome\nand another\nsubject's diet")) +
  scale_fill_manual(values = c("#5f86b7","#5a2071","#64baaa")) +
  geom_errorbar(data = summary.master, aes(x = variable, ymin = mean - se, ymax = mean + se), width = 0.4, inherit.aes = F) +
  theme_classic() +
  ylab("Prediction accuracy\n(Pearson correlation of species abundance)") +
  xlab("Model training set") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7)) 
  #geom_text(data = summary.master, aes(x = variable, y = mean + se, label=groups), vjust=-1, size = 3)

# pring the plot
mybarplot

ggsave("figures/prediction_bar.pdf", mybarplot, height = 2.5, width = 2.5, device = "pdf")

  # make paired difference plots by person

require(gridExtra)
require(directlabels)

mb.dt <- master.df.melt[master.df.melt$variable %in% c("own.mb.corr", "own.diet.corr"),]

mb.dt.wide <- merge(master.df, own.diet.pvalues.melt)
mb.dt.wide <- merge(mb.dt.wide, own.diet.pc1.pvals.melt)
mb.dt.wide <- merge(mb.dt.wide, own.diet.pc2.pvals.melt)
mb.dt.wide <- merge(mb.dt.wide, own.diet.pc3.pvals.melt)

mb.dt.wide <- mb.dt.wide[,!colnames(mb.dt.wide) %in% c("others.diet.corr")]

mb.dt.wide$sig<- ifelse(mb.dt.wide$own.diet.p <= 0.25,"sig", "ns")

mb.dt.wide$groups <- ifelse(mb.dt.wide$sig == "sig", mb.dt.wide$Taxa, NA)

mb.dt.for.merge <- mb.dt.wide[colnames(mb.dt.wide) %in% c("Taxa", "UserName", "groups")]

mb.dt <- merge(mb.dt, mb.dt.for.merge)


# make multiple plots in a loop

single_plots <- as.list(NULL)

plot <- mb.dt[order(mb.dt$UserName, decreasing = F),]
plot <- plot[order(plot$groups, na.last = F),]

group.names <- unique(mb.dt$groups)
group.colors <- c("#40655e", "#bd836e", "#48bea1", "#972554", "#79ac3d", "#6029ce", "#fb5de7", "#9ca6f4", "#9f2108", 
                  "#30b6f3", "#5b3c5e", "#dd8eeb", "#31c52f", "#fa5a45", "#9cb0a2", "#1f3ca6", "#ef972d", "#056e12", 
                  "#7d4400","#858c69","#008c4b", "#00138c","#5b430b", "#8c7723", "#468c75","#8091ff", "#ff80c4","red", 
                  "#5940ff","#b70d61","#462881", "#a50fa9", "#167b2b", "#fe8f06","maroon","darkgrey", "blue", "black")


names(group.colors) <- group.names



# view colors
cols <- function(a) image(1:31, 1, as.matrix(1:31), col=a, axes=FALSE , xlab="", ylab="")
cols(group.colors)

for (i in 1:length(unique(plot$UserName))) {
  j <- as.character(unique(plot$UserName)[[i]])
  plot_un <- plot[plot$UserName == j,]
  
  plot_un <- plot_un[order(plot_un$groups, na.last = T),]
  
  single_plots[[i]] <- ggplot(data = plot_un, aes(x = variable, y = value)) + 
    facet_free(.~UserName, scales = "free") +
    geom_line(aes(group = Taxa), color = "grey", alpha = 0.5) + 
    geom_line(aes(group = Taxa, color = factor(groups)), size = 1) +
    scale_color_manual(values = group.colors) +
    #geom_dl(aes(label = groups), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.5, hjust = 1)) +
    geom_dl(aes(label = groups, color = groups), method = "last.polygons") +
    geom_point() +
    theme_classic() +
    ylim(0,1) +
    scale_x_discrete(expand = c(0,0.1)) +
    expand_limits(x = 2.5) +
    theme(legend.position = "none",
          legend.text = element_text(size = 7),
          legend.title = element_blank(),
          axis.title = element_blank()) 
}


# get find out what people have significant diet-related species
keep_plots <- unique(subset(plot, !is.na(plot$groups))$UserName)
#keep_plots <- keep_plots[!keep_plots %in% c("11", "12")] # don't plot soylent people
loc <- which(mynames %in% keep_plots)


# plot just those graphs
grid.arrange(grobs = single_plots[loc], nrow = 1) # nrow=round(length(loc))
g <-arrangeGrob(grobs = single_plots[loc], nrow=1)  


ggsave("figures/species.person.model.png", g, device = "png", height = 3, width = 9)


#### make plots of pcs and pvals by species
require(ggbeeswarm)
modelpcs <- merge(own.diet.pc1.melt, own.diet.pc1.pvals.melt)
modelpcs <- merge(modelpcs, own.diet.pc2.melt)
modelpcs <- merge(modelpcs, own.diet.pc2.pvals.melt)
modelpcs <- merge(modelpcs, own.diet.pc3.melt)
modelpcs <- merge(modelpcs, own.diet.pc3.pvals.melt)

modelpcs$pc1sig <- ifelse(modelpcs$own.diet.pc1.pvals < 0.25, "sig", "ns")
modelpcs$pc2sig <- ifelse(modelpcs$own.diet.pc2.pvals < 0.25, "sig", "ns")
modelpcs$pc3sig <- ifelse(modelpcs$own.diet.pc3.pvals < 0.25, "sig", "ns")

# reorder taxa for plotting

load(file = "data/test_personal_diet_taxastrings.Rdata")
taxasplit <- colsplit(taxonomy$taxastring, ";", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
taxonomy <- cbind(taxonomy, taxasplit)
taxonomy$Phylum <- ifelse(is.na(taxonomy$Phylum), taxonomy$Kingdom, taxonomy$Phylum)
taxonomy$Class <- ifelse(is.na(taxonomy$Class), taxonomy$Phylum, taxonomy$Class)
taxonomy$Order <- ifelse(is.na(taxonomy$Order), taxonomy$Class, taxonomy$Order)
taxonomy$Family <- ifelse(is.na(taxonomy$Family), taxonomy$Order, taxonomy$Family)
taxonomy$Genus <- ifelse(is.na(taxonomy$Genus), taxonomy$Family, taxonomy$Genus)


modelpcs <- merge(modelpcs, taxonomy, by.x = "Taxa", by.y = "mbname")

# reorder
test <- modelpcs[with(modelpcs, order(Phylum, Class, Order, Family, Genus)),]
taxaorder <- unique(test$Taxa)
familyorder <- unique(test$Family)
phylumorder <- unique(test$Phylum)

modelpcs$Taxa <- factor(modelpcs$Taxa, levels = taxaorder, ordered = T)
modelpcs$Family <- factor(modelpcs$Family, levels = familyorder, ordered = T)
modelpcs$Phylum <- factor(modelpcs$Phylum, levels = phylumorder, ordered = T)

pcplot <- ggplot(data = modelpcs) +
  geom_quasirandom(data = subset(modelpcs,pc1sig ="ns"), aes(x = Taxa, y = own.diet.pc1), size = 2, alpha = 0.5, color = "grey")+
  geom_point(data = subset(modelpcs, pc1sig == "sig"), aes(x = Taxa, y = own.diet.pc1), size = 2, alpha =1, color = "#5a2071") +
  #geom_label(data = subset(modelpcs, pc2sig == "sig"), aes(label = UserName), size = 2, nudge_y = 7, nudge_x = 1) +
  scale_y_continuous(limits = c(-30,30), breaks = c(-20,-10,0,10,20)) + 
  theme_classic() +
  theme(#axis.text.x = element_text(angle =60, hjust = 1),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  facet_grid(.~Phylum + Family + Genus, space = "free", scale = "free") +
  ylab("Own Diet PC1") +
  xlab("")


pcplot2 <- ggplot(data = modelpcs) +
  geom_quasirandom(data = subset(modelpcs,pc2sig ="ns"), aes(x = Taxa, y = own.diet.pc2), size = 2, alpha = 0.5, color = "grey")+
  geom_point(data = subset(modelpcs, pc2sig == "sig"), aes(x = Taxa, y = own.diet.pc2), size = 2, alpha =1, color = "#5a2071") +
  #geom_label(data = subset(modelpcs, pc2sig == "sig"), aes(label = UserName), size = 2, nudge_y = 7, nudge_x = 1) +
  scale_y_continuous(limits = c(-30,30), breaks = c(-20,-10,0,10,20)) + 
  theme_classic() +
  theme(#axis.text.x = element_text(angle =60, hjust = 1),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  facet_grid(.~Phylum + Family + Genus, space = "free", scale = "free") +
  ylab("Own Diet PC2") +
  xlab("")


pcplot3 <- ggplot(data = modelpcs) +
  geom_quasirandom(data = subset(modelpcs,pc3sig ="ns"), aes(x = Taxa, y = own.diet.pc3), size = 2, alpha = 0.5, color = "grey")+
  geom_point(data = subset(modelpcs, pc3sig == "sig"), aes(x = Taxa, y = own.diet.pc3), size = 2, alpha =1, color = "#5a2071") +
  #geom_label(data = subset(modelpcs, pc3sig == "sig"), aes(label = UserName), size = 2, nudge_y = 7, nudge_x = 1) +
  scale_y_continuous(limits = c(-30,30), breaks = c(-20,-10,0,10,20)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle =60, hjust = 1),
        panel.grid.major = element_line(colour = "grey"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  facet_grid(.~Phylum + Family + Genus, space = "free", scale = "free") +
  ylab("Own Diet PC3") +
  xlab("")



ggsave("figures/pc1.png",plot = pcplot,width = 11, height = 2.7, device = "png")
ggsave("figures/pc2.png", plot = pcplot2, width = 11, height = 2.7, device = "png")
ggsave("figures/pc3.png", plot = pcplot3, width = 11, height = 4.2,device = "png")




# check if any speices are significantly influenced by diet 
test1 <- apply(own.diet.pc1[,2:ncol(own.diet.pc1)], 1, function(x) t.test(x)$p.value)
test1 <- p.adjust(test1, "fdr")
sum(test1 < 0.25)

test2 <- apply(own.diet.pc2[,2:ncol(own.diet.pc2)], 1, function(x) t.test(x)$p.value)
test2 <- p.adjust(test2, "fdr")
sum(test2 < 0.25)

test3 <- apply(own.diet.pc3[,2:ncol(own.diet.pc3)], 1, function(x) t.test(x)$p.value)
test3 <- p.adjust(test3, "fdr")
sum(test3 < 0.25)




