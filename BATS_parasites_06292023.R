## Code for processing 18S rRNA gene metabarcoding data (V4) from BATS 
# Updated 11/8/2023
# Sean R. Anderson

# Load packages
library(tidyverse)
library(vegan)
library(qiime2R)
library(phyloseq)
library(fantaxtic)
library(RColorBrewer)
library(microbiome)
library(factoextra)
library(Matrix)
library(ggpubr)
library(SpiecEasi)
library(microeco)
library(file2meco)
library(mdatools)
library(ComplexUpset)
library(reshape)
library(reshape2)
library(pairwiseAdonis)
library(magrittr)

# Load 18S taxonomy file - assigned with new PR2 database release (V 5.0.1)
taxonomy <- read_qza(file="18S-tax_new.qza")
tax_tab <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)

# Load ASV count table 
table <- read_qza(file = "table.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 

# Load 18S metadata file from BATS
metadata <- read.table("metadata_v2.txt", header=TRUE, row.names=1, check.names=F, sep="\t") 

# Merge into phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax_tab)), otu_table(count_tab, taxa_are_rows = T), sample_data(metadata))

# Remove unwanted groups
ps_new = subset_taxa(ps, Division !="Streptophyta" |is.na(Division))
ps_new = subset_taxa(ps_new, Division !="Rhodophyta" |is.na(Division))
ps_new = subset_taxa(ps_new, Domain !="Bacteria" |is.na(Domain))
ps_new <- subset_taxa(ps_new, Subdivision!="Opisthokonta_X"|is.na(Subdivision))
ps_new <- subset_taxa(ps_new, Class!="Unassigned", Prune = T)
ps_new = subset_taxa(ps_new, Class!="Craniata"|is.na(Class))

ps_new = name_na_taxa(ps_new) # Adds an unassigned label to better identify lowest possible taxonomic assignment

# Remove samples with <5,000 reads and unwanted depths (100 and 400 m)
ps_sub = prune_samples(sample_sums(ps_new)>=5000, ps_new)
ps_sub <- subset_samples(ps_sub, Nominal_Depth != "100" & Nominal_Depth != "400")

# Remove singletons (ASVs present only once)
ps_filt = filter_taxa(ps_sub, function (x) {sum(x) > 1}, prune=TRUE)

# Rarefaction curves for 18S - Figure SI 10
rare_18S <- ggrare(ps_filt, step = 100, plot = FALSE, parallel = FALSE, se = FALSE)
rare_18S + theme(legend.position = "none") + theme_bw()+ theme(legend.position = "right") + 
  facet_wrap(~Nominal_Depth, scales="free_y",nrow=4,ncol=3)
ggsave(filename = "18S_rare.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 7, height = 5, dpi = 150)

# Estimate minimum, mean, and maximum 18S read counts
ps_min <- min(sample_sums(ps_filt))
ps_mean <- mean(sample_sums(ps_filt))
ps_max <- max(sample_sums(ps_filt))

# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps_filt, sample.size = min(sample_sums(ps_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Rename ASVs in sequential order
taxa_names(ps_rare) <- paste0("ASV", seq(ntaxa(ps_rare)))

# Export 18S ASV information - filtered, rarefied, and relabeled ASVs (used in data analysis for this study) - Table S4
OTU_filt = as(otu_table(ps_rare), "matrix")
TAX_filt = as(tax_table(ps_rare), "matrix")
merge_18S_filt <- cbind(OTU_filt,TAX_filt)
write.csv(merge_18S_filt, file="Table_S4.csv", row.names=T)

# Plot relative abundance of 18S groups with depth at BATS 
dataset <- phyloseq2meco(ps_rare) # Convert to microeco object
dataset1 <- dataset$merge_samples(use_group = "Nominal_Depth") # Group the data by depth
t1 <- trans_abund$new(dataset = dataset1, taxrank = "Class", ntaxa = 10) # Create an abundance object with the top 10 class level 18S groups

# Plot 18S class level relative abundance - Figure 1A
g1 <- t1$plot_bar(bar_type = "full", use_alluvium = FALSE, clustering = FALSE,barwidth = NULL, xtext_size = 12, color_values = c("#332288","#CC6677","#DDCC77","#117733","#88CCEE", "#882255", "#44AA99", "#999933","#AA4499" , "#EB7357"),others_color = "#757575",
                  order_x = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))
g1 + coord_flip() + geom_col(colour = "black") + labs(x="Depth (m)",y="Relative abundance (%)")+theme(text = element_text(size = 12))
ggsave(filename = "syn_stacked_class.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 5, dpi = 300)

# Relative abundance plots within Syndiniales (order and family levels) 
ps_Syn =subset_taxa(ps_rare, Class=="Syndiniales",Prune = T) # Subset to only Syndiniales
dataset <- phyloseq2meco(ps_Syn) 
dataset1 <- dataset$merge_samples(use_group = "Nominal_Depth")
t1 <- trans_abund$new(dataset = dataset1, taxrank = "Order", ntaxa = 10) # Create an abundance object with the top 10 order level groups within Syndiniales

# Plot order level Syndiniales relative abundance - Figure SI 2
g1 <- t1$plot_bar(bar_type = "full", use_alluvium = FALSE, clustering = FALSE,barwidth = NULL, xtext_size = 12, color_values = c("#A5CFCC","#0E899F","#F1B6A1", "#D4A52A", "#E3E5DB", "#A83860", "#ED91BC"),others_color = "#757575",
                  order_x = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))
g1 + coord_flip()+geom_col(colour = "black")+ labs(x="Depth (m)",y="Relative abundance (%)")+theme(text = element_text(size = 12))
ggsave(filename = "syn_stacked_order.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 4.5, height = 5, dpi = 300)

# Plot family level Syndiniales relative abundance - Figure 1B
t1 <- trans_abund$new(dataset = dataset1, taxrank = "Family", ntaxa = 10) 
g1 <- t1$plot_bar(bar_type = "full", use_alluvium =FALSE, clustering = FALSE,barwidth = NULL, xtext_size = 12, color_values = c("#4477AA","#EE6677","#228833", "#CCBB44", "#66CCEE", "#AA3377", "#EE7733","#CC3311","#009988","#EE3377"),others_color = "#757575",
                  order_x = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))
g1 + coord_flip()+geom_col(colour = "black")+ labs(x="Depth (m)",y="Relative abundance (%)")+theme(text = element_text(size = 12))
ggsave(filename = "syn_stacked_family.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5.5, height = 5, dpi = 300)

# Plot 18S stacked bar plots with depth and across seasons - Figure SI 1
barplot <- ps_rare %>%
  tax_glom(taxrank = "Class", NArm=FALSE) %>% # Agglomerate at class level
  transform_sample_counts(function(OTU) 100* OTU/sum(OTU)) %>% # Transform to relative abundance
  psmelt() %>% # Melt data
  group_by(Nominal_Depth,Class,Season) %>% # Group by depth and season
  summarise_at("Abundance", .funs = mean) # Summarize at the mean

barplot_18S = barplot # Rename
barplot_18S$Class<- as.character(barplot_18S$Class)
barplot_18S$Class[barplot_18S$Abundance < 2]= "Other" # Set a threshold where groups that are < 2% relative abundance are grouped into "Other" category
barplot_18S$Nominal_Depth = as.factor(barplot_18S$Nominal_Depth)
barplot_18S[barplot_18S== "Annelida" ] <- "Other" # Remove other less abundant groups
barplot_18S[barplot_18S==  "Mamiellophyceae"] <- "Other"
barplot_18S[barplot_18S==  "Pelagophyceae"] <- "Other"
barplot_18S[barplot_18S==  "Urochordata"] <- "Other"
p <- ggplot(data=barplot_18S, aes(x=fct_rev(Nominal_Depth), y=Abundance, fill=Class))
p$data$Class <- factor(p$data$Class, levels = c("Other","Prymnesiophyceae","Sagenista","RAD-A", "Cnidaria","RAD-B","Acantharea","Arthropoda", "Polycystinea", "Dinophyceae", "Syndiniales")) # Set order of class groups in the plot
p$data$Season <- factor(p$data$Season, levels = c("Winter","Spring","Summer","Fall")) # Order the seasons
p + geom_bar(aes(), stat="identity", position="fill", width = 0.9)+
  scale_y_continuous(expand = c(0, 0))+
  geom_hline(yintercept=0) + theme_bw()+ 
  scale_fill_manual(values= rev(c("#332288","#CC6677","#DDCC77","#117733","#88CCEE", "#882255", "#44AA99", "#999933","#AA4499" , "#EB7357","#757575")))+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + 
  theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=11, ncol=1)) + 
  facet_wrap(~Season,ncol=2,nrow=2,)  + 
  coord_flip() +labs(y = "Relative abundance (%)", x = "Depth (m)") 
ggsave(filename = "18S_barplots_time.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

# Beta diversity for all Syn samples in the dataset - Figure 1C
ordu = ordinate(ps_Syn, "PCoA", "bray") # PCoA ordination based on Bray-Curtis
p = plot_ordination(ps_Syn, ordu, color="Nominal_Depth")
p$data$Nominal_Depth <- as.factor(p$data$Nominal_Depth) # Convert sampling depth to factor
mycolors = c(brewer.pal(name="Paired", n = 12)) # Set color palette for the plot
p +theme_bw() + scale_fill_manual(values=mycolors) +
  geom_point(aes(fill=Nominal_Depth),size = 6, shape = 21, colour = "black")  +   
  theme(text = element_text(size=14)) 
ggsave(filename = "PCoA_syn.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, dpi = 300)

# PERMANOVA to test significance of sampling depth (or season)
metadata <- as(sample_data(ps_Syn), "data.frame")
metadata$Nominal_Depth=as.factor(metadata$Nominal_Depth)
metadata$Season=as.factor(metadata$Season)
bray = phyloseq::distance(ps_Syn, method="bray") 
adonis2(phyloseq::distance(ps_Syn, method="bray")~Nominal_Depth, data = metadata,perm=9999)
pairwise.adonis(as.dist(bray), as.factor(metadata$Nominal_Depth), p.adjust.m = 'bonferroni',perm = 9999) # Include in Table S1

# Run dispersion tests on sampling depth (or season)
perm.eg.betadisper <- betadisper(bray, group = metadata$Nominal_Depth, type = "centroid")
anova(perm.eg.betadisper)
adonis2(dist(perm.eg.betadisper$distances) ~ metadata$Nominal_Depth,perm=9999)

# Beta diversity variation with depth - Figure 1D
dataset <- phyloseq2meco(ps_Syn)
dataset$cal_betadiv(unifrac = FALSE)
t1 <- trans_beta$new(dataset = dataset, group = "Nominal_Depth", measure = "bray")
t1$cal_group_distance(within_group = T)
t1$cal_group_distance_diff(method = "anova") # Significance with ANOVA
g1 <- t1$plot_group_distance(boxplot_add = "mean",color_values = "#DDAA33",xtext_keep = TRUE)
g1$data$Nominal_Depth <- factor(g1$data$Nominal_Depth, levels = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))
g1 + geom_boxplot(fill="#DDAA33")+scale_x_discrete(limits = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1, size=12)) +coord_flip()
ggsave(filename = "Bray_Syn.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 4, height = 5, dpi = 300)

# Alpha diversity between photic and aphotic zones - Figure 1E
ps_Syn2 <- subset_taxa(ps_sub, Class=="Syndiniales",Prune = T) # Subset to only Syndiniales - include singletons for diversity analysis
rich_syn <- estimate_richness(ps_Syn2, measures=c("Observed", "Shannon")) # Estimate richness and Shannon diversity index
zone = sample_data(ps_Syn2)$vert_bin # Define vertical zones (photic vs. aphotic)
rich_all <- data.frame(rich_syn,zone) # Merge diversity values with vertical zones
rich_all$zone = as.factor(rich_all$zone)
p <- ggplot(rich_all, aes(x=zone, y=Shannon, fill=zone))
p$data$zone <- factor(p$data$zone, levels = c("Photic","Aphotic")) # Order depth zones
p + geom_boxplot(outlier.shape = NA) + theme_bw() + 
  theme(text = element_text(size=14)) + ylab("Syndiniales shannon index") + theme(legend.position="right")+ scale_fill_manual(values=c("#DDAA33","#004488"))+
  geom_point(aes(fill=zone), size =5, shape = 21, colour = "black", position=position_jitterdodge(0.4))+
  theme(axis.text.x=element_text(angle=45, hjust=1, color="black"))+
  theme(axis.title.x =element_blank()) + 
  stat_compare_means(method= "wilcox.test",comparisons = list(c("Photic", "Aphotic")), label = "p.signif", exact = TRUE) # Compare means with Wilcoxon test
ggsave(filename = "Shan_Syn.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 5, dpi = 300)

### PLSR analysis with Syndiniales ASVs over all depths at BATS
class_18S <- tax_glom(ps_rare, taxrank = "Class",NArm=FALSE) # Agglomerate to class level
x1 = psmelt(class_18S) # Melt the data

# Group the data based on read counts at the class level and across environmental variables
mean.18S = x1 %>% 
  dplyr::group_by(Sample,Class,Nominal_Depth) %>%
  dplyr::summarise_at(.vars = c("Temp", "CTD_S", "oxygen", "PO4", "NO3", "NO2", "SiO4", "POC","PON","CO2","Alk", "TOC", "TN", "Bact", "Abundance"), .funs = mean)

# Subset to Syndiniales ASVs
abund_18S <- subset(mean.18S, grepl("Syndiniales", mean.18S$Class)) 

# Prepare for PLSR (split environmental factors from read counts)
myData = na.omit(abund_18S) # Remove missing variables
factors = myData[, c(1,4:17)] # Split environmental variables
factors$Sample <- NULL # Remove sample column
abund_split = myData[, c(2,18)] # Split Syndiniales read counts
abund_final = abund_split[, c(2)] # Get rid of the first column

# Run PLSR with Syndiniales read counts
m1 = pls(factors, abund_final,7, scale = TRUE,center = TRUE, cv = 1, ncomp.selcrit = "min") # Run the core pls function
show(m1$ncomp.selected) # 4 components are optimal 
plotRMSE(m1) # Root mean squared plot
summary(m1) # Summary information for the model
summary(m1$res$cal) # How much variation is explained by components
plotPredictions(m1) # Plot predictions
plotPredictions(m1$res$cal, ncomp = 4, show.stat = TRUE)
plotVIPScores(m1) # Plot VIP scores for the different environmental factors
plotVIPScores(m1, ncomp = 4, type = "h", show.labels = TRUE)

# Remove outliers in the dataset and re-run the PLSR
m1 = pls(factors, abund_final,7, scale = TRUE, center=TRUE, cv = 1, lim.type = "ddmoments")
plotXYResiduals(m1, show.labels = TRUE, labels = "indices")
m = setDistanceLimits(m1, lim.type = "ddrobust") # Use robust method to identify outliers
plotXYResiduals(m, show.labels = FALSE, labels = "indices") # Plot outliers - Figure SI 5

outliers = which(categorize(m, m$res$cal) == "outlier") # Get row indices for outliers in calibration set
Xo = factors[outliers, , drop = FALSE] # Keep data for outliers in separate matrices
yo = factors[outliers, , drop = FALSE]
X = factors[-outliers, , drop = FALSE] # Remove the rows with outliers from the data
y = factors[-outliers, , drop = FALSE]

X0 = abund_final[outliers, , drop = FALSE] 
y0 = abund_final[outliers, , drop = FALSE]
x = abund_final[-outliers, , drop = FALSE]
Y = abund_final[-outliers, , drop = FALSE]

m1 = pls(X, x,7, scale = TRUE, center=TRUE, cv = 1, ncomp.selcrit = "min") # Re-run the PLSR without outliers
plotRMSE(m1) # Plot RMSE - Figure SI 5
show(m1$ncomp.selected) # Show number of optimal components
plotPredictions(m1) # Plot predictions
plotPredictions(m1$res$cal, ncomp = 4, show.stat = TRUE) # Plot predictions with optimal comp - Figure SI 5
plotPredictions(m1, ncomp = 4)
vip2 = vipscores(m1, ncomp = 4) # Save the VIP values
summary(m1) # Summary of the new model
plot(m1$coeffs, ncomp = 4,  show.labels = TRUE) # Plot model coefficients - Figure SI 5

# Plot the VIP values from PLSR model - Figure 1F
VIP_all = as.data.frame(vip2) # Convert to data frame 
names(VIP_all)[1] <- "Syndiniales" # Relabel column
VIP_all$VAR <- row.names(VIP_all) # Create new column with environmental variables
p <- ggplot(VIP_all, aes(x=reorder(VAR, -Syndiniales), y=Syndiniales))
p + geom_col(show.legend = TRUE,position = "dodge",colour="black",fill="darkgray") +
  theme_bw()+ geom_hline(yintercept=1, linetype=2, colour="black")+
  tidytext::scale_x_reordered() + labs(y="VIP Scores")+theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0),limits=c(0,3)) + theme(axis.text.x=element_text(angle=45, hjust=1, color="black"))
ggsave(filename = "PLS_Syndiniales.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 4, dpi = 300)

# Plot oxygen concentration vs. Syndiniales relative abundance across all samples 
class_18S_rel <-  transform_sample_counts(class_18S, function(x) x / sum(x)) # Convert to relative abundance
x1 = psmelt(class_18S_rel) # Melt data
x1$Abundance <- x1$Abundance * 100 # Multiply to get percent

# Repeat grouping by factors and subset to Syndiniales
mean.18S = x1 %>% 
  dplyr::group_by(Sample,Class,Nominal_Depth) %>%
  dplyr::summarise_at(.vars = c("Temp", "CTD_S", "oxygen", "PO4", "NO3", "NO2", "SiO4", "POC","PON","CO2","Alk", "TOC", "TN", "Bact", "Abundance"), .funs = mean)

abund_18S <- subset(mean.18S, grepl("Syndiniales", mean.18S$Class)) # Subset to Syndiniales - based on relative abundance now
myData = na.omit(abund_18S) # Remove missing values
myData$Nominal_Depth = as.factor(myData$Nominal_Depth) # Convert depth to factor

# Oxygen plotted vs. Syndiniales - can substitute other factors - Figure 1G
p <- ggscatter(myData, x = "oxygen", y = "Abundance", cor.method = "spearman", cor.coef =T, size = 3, add.params = list(color = "black", fill = "darkgray",col="Nominal_Depth"), add = "reg.line", conf.int = TRUE, ylab = "Syndiniales relative abundance (%)", xlab = "oxygen (µmol/kg)") +
  geom_point(aes(fill=Nominal_Depth),size = 5, shape = 21, colour = "black") +scale_fill_manual(values=mycolors) + theme(legend.position = "right") + theme(text = element_text(size = 18))  + theme_bw()           
p
ggsave(filename = "Syndiniales_oxygen.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 4, dpi = 300)

# Plot environmental variables with sampling depth at BATS  - Figure SI 4
meta1 = myData[,c(3:17)] # Select columns of interest
new = reshape2::melt(meta1,id='Nominal_Depth')
new$Nominal_Depth = as.factor(new$Nominal_Depth)
p <- ggplot(new, aes(x=fct_rev(Nominal_Depth), y=value))
p$data$variable <- factor(p$data$variable, levels =c("Temp","CTD_S","oxygen","PO4", "NO3","NO2","SiO4","CO2","Alk","TN","POC","PON","TOC","Bact"))
p + theme_bw() + geom_boxplot()+stat_boxplot(geom = "errorbar", width = 0.5)+
  theme(text = element_text(size=14)) +theme(legend.position="right")+ 
  labs(x="Depth (m)")+ theme(axis.text.x=element_text(angle=45, hjust=1, color="black"))+coord_flip()+ facet_wrap(~variable,scales="free_x",ncol=4)
ggsave(filename = "factors_depth_new.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 5, height = 8, dpi = 300)

### Prepare for network analysis in SPIEC-EASI - partitioning data into twelve depths
ps_18S_surface <- subset_samples(ps_rare, Nominal_Depth == "1")
ps_18S_40 <- subset_samples(ps_rare, Nominal_Depth == "40")
ps_18S_80 <- subset_samples(ps_rare, Nominal_Depth == "80")
ps_18S_120 <- subset_samples(ps_rare, Nominal_Depth == "120")
ps_18S_160 <- subset_samples(ps_rare, Nominal_Depth == "160")
ps_18S_200 <- subset_samples(ps_rare, Nominal_Depth == "200")
ps_18S_250 <- subset_samples(ps_rare, Nominal_Depth == "250")
ps_18S_300 <- subset_samples(ps_rare, Nominal_Depth == "300")
ps_18S_500 <- subset_samples(ps_rare, Nominal_Depth == "500")
ps_18S_600 <- subset_samples(ps_rare, Nominal_Depth == "600")
ps_18S_800 <- subset_samples(ps_rare, Nominal_Depth == "800")
ps_18S_1000 <- subset_samples(ps_rare, Nominal_Depth == "1000")

# Subset to the top 150 most abundant ASVs in each depth
filter_18S_surface = prune_taxa(names(sort(taxa_sums(ps_18S_surface), TRUE))[1:150], ps_18S_surface) 
filter_18S_40 = prune_taxa(names(sort(taxa_sums(ps_18S_40), TRUE))[1:150], ps_18S_40)
filter_18S_80 = prune_taxa(names(sort(taxa_sums(ps_18S_80), TRUE))[1:150], ps_18S_80)
filter_18S_120 = prune_taxa(names(sort(taxa_sums(ps_18S_120), TRUE))[1:150], ps_18S_120)
filter_18S_160 = prune_taxa(names(sort(taxa_sums(ps_18S_160), TRUE))[1:150], ps_18S_160)
filter_18S_200 = prune_taxa(names(sort(taxa_sums(ps_18S_200), TRUE))[1:150], ps_18S_200)
filter_18S_250 = prune_taxa(names(sort(taxa_sums(ps_18S_250), TRUE))[1:150], ps_18S_250)
filter_18S_300 = prune_taxa(names(sort(taxa_sums(ps_18S_300), TRUE))[1:150], ps_18S_300)
filter_18S_500 = prune_taxa(names(sort(taxa_sums(ps_18S_500), TRUE))[1:150], ps_18S_500)
filter_18S_600 = prune_taxa(names(sort(taxa_sums(ps_18S_600), TRUE))[1:150], ps_18S_600)
filter_18S_800 = prune_taxa(names(sort(taxa_sums(ps_18S_800), TRUE))[1:150], ps_18S_800)
filter_18S_1000 = prune_taxa(names(sort(taxa_sums(ps_18S_1000), TRUE))[1:150], ps_18S_1000)

# Format OTU tables for analysis
surface.otu <- t(data.frame(phyloseq::otu_table(filter_18S_surface), check.names = FALSE)) 
otu.40 <- t(data.frame(phyloseq::otu_table(filter_18S_40), check.names = FALSE)) 
otu.80 <- t(data.frame(phyloseq::otu_table(filter_18S_80), check.names = FALSE)) 
otu.120 <- t(data.frame(phyloseq::otu_table(filter_18S_120), check.names = FALSE)) 
otu.160 <- t(data.frame(phyloseq::otu_table(filter_18S_160), check.names = FALSE)) 
otu.200 <- t(data.frame(phyloseq::otu_table(filter_18S_200), check.names = FALSE)) 
otu.250 <- t(data.frame(phyloseq::otu_table(filter_18S_250), check.names = FALSE)) 
otu.300 <- t(data.frame(phyloseq::otu_table(filter_18S_300), check.names = FALSE)) 
otu.500 <- t(data.frame(phyloseq::otu_table(filter_18S_500), check.names = FALSE)) 
otu.600 <- t(data.frame(phyloseq::otu_table(filter_18S_600), check.names = FALSE)) 
otu.800 <- t(data.frame(phyloseq::otu_table(filter_18S_800), check.names = FALSE)) 
otu.1000 <- t(data.frame(phyloseq::otu_table(filter_18S_1000), check.names = FALSE)) 

# Run each depth network using the same parameters 
se <- spiec.easi(surface.otu, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.40 <- spiec.easi(otu.40, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.80 <- spiec.easi(otu.80, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.120 <- spiec.easi(otu.120, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.160 <- spiec.easi(otu.160, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.200 <- spiec.easi(otu.200, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.250 <- spiec.easi(otu.250, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.300 <- spiec.easi(otu.300, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.500 <- spiec.easi(otu.500, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.600 <- spiec.easi(otu.600, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.800 <- spiec.easi(otu.800, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 
se.1000 <- spiec.easi(otu.1000, method ='mb', nlambda = 40, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05)) 

# Function to get weights and edge data into correct format for exporting. Re-run this function for each network output above. 
convertSEToTable <- function(se.1000,sp.names){ # Run this function to get the edge weights
  sebeta <- symBeta(getOptBeta(se.1000), mode='maxabs') 
  elist     <- summary(sebeta)
  elist[,1] <- sp.names[elist[,1]]
  elist[,2] <- sp.names[elist[,2]]
  elist[,4] <- paste(elist[,1],elist[,2])
  full_e <- expand.grid(sp.names,sp.names)
  rownames(full_e) <- paste(full_e[,1],full_e[,2])
  full_e[,"Weight"] <- 0
  full_e[elist[,4],"Weight"] <- elist[,3]
  x <- expand.grid(1:length(sp.names),1:length(sp.names))
  full_e[x[,"Var1"]>x[,"Var2"],"Weight"] <- NA
  return(as.data.frame(full_e,stringsAsFactors=F))
}

tab.se <- convertSEToTable(se.1000, sp.names=colnames(otu.1000)) 
tab.se.filtered <- tab.se %>% filter(is.finite(Weight))

# Remove values that are 0 for each network
filt_1 <- subset(tab.se.filtered, abs(Weight) > 0) 
filt_40 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_80 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_120 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_160 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_200 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_250 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_300 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_500 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_600 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_800 <- subset(tab.se.filtered, abs(Weight) > 0)
filt_1000 <- subset(tab.se.filtered, abs(Weight) > 0)

# Merge network and taxa files for each depth network and then merge all together into giant table
taxALL_1 <- tax_table(filter_18S_surface)
taxALL_1 <- as.data.frame(taxALL_1, check.names=FALSE)
taxALL_1 <- cbind(rownames(taxALL_1), data.frame(taxALL_1, row.names=NULL))
colnames(taxALL_1)[1] <- "Var1"
newtab_1 <- merge(filt_1, taxALL_1[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_1)[1] <- "Var2"
newtab_final_1 <- merge(newtab_1, taxALL_1[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_1 = newtab_final_1 %>%
  mutate(Network = "1m")

taxALL_40 <- tax_table(filter_18S_40)
taxALL_40 <- as.data.frame(taxALL_40, check.names=FALSE)
taxALL_40 <- cbind(rownames(taxALL_40), data.frame(taxALL_40, row.names=NULL))
colnames(taxALL_40)[1] <- "Var1"
newtab_40 <- merge(filt_40, taxALL_40[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_40)[1] <- "Var2"
newtab_final_40 <- merge(newtab_40, taxALL_40[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_40 = newtab_final_40 %>%
  mutate(Network = "40m")

taxALL_80 <- tax_table(filter_18S_80)
taxALL_80 <- as.data.frame(taxALL_80, check.names=FALSE)
taxALL_80 <- cbind(rownames(taxALL_80), data.frame(taxALL_80, row.names=NULL))
colnames(taxALL_80)[1] <- "Var1"
newtab_80 <- merge(filt_80, taxALL_80[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_80)[1] <- "Var2"
newtab_final_80 <- merge(newtab_80, taxALL_80[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_80 = newtab_final_80 %>%
  mutate(Network = "80m")

taxALL_120 <- tax_table(filter_18S_120)
taxALL_120 <- as.data.frame(taxALL_120, check.names=FALSE)
taxALL_120 <- cbind(rownames(taxALL_120), data.frame(taxALL_120, row.names=NULL))
colnames(taxALL_120)[1] <- "Var1"
newtab_120 <- merge(filt_120, taxALL_120[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_120)[1] <- "Var2"
newtab_final_120 <- merge(newtab_120, taxALL_120[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_120 = newtab_final_120 %>%
  mutate(Network = "120m")

taxALL_160 <- tax_table(filter_18S_160)
taxALL_160 <- as.data.frame(taxALL_160, check.names=FALSE)
taxALL_160 <- cbind(rownames(taxALL_160), data.frame(taxALL_160, row.names=NULL))
colnames(taxALL_160)[1] <- "Var1"
newtab_160 <- merge(filt_160, taxALL_160[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_160)[1] <- "Var2"
newtab_final_160 <- merge(newtab_160, taxALL_160[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_160 = newtab_final_160 %>%
  mutate(Network = "160m")

taxALL_200 <- tax_table(filter_18S_200)
taxALL_200 <- as.data.frame(taxALL_200, check.names=FALSE)
taxALL_200 <- cbind(rownames(taxALL_200), data.frame(taxALL_200, row.names=NULL))
colnames(taxALL_200)[1] <- "Var1"
newtab_200 <- merge(filt_200, taxALL_200[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_200)[1] <- "Var2"
newtab_final_200 <- merge(newtab_200, taxALL_200[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_200 = newtab_final_200 %>%
  mutate(Network = "200m")

taxALL_250 <- tax_table(filter_18S_250)
taxALL_250 <- as.data.frame(taxALL_250, check.names=FALSE)
taxALL_250 <- cbind(rownames(taxALL_250), data.frame(taxALL_250, row.names=NULL))
colnames(taxALL_250)[1] <- "Var1"
newtab_250 <- merge(filt_250, taxALL_250[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_250)[1] <- "Var2"
newtab_final_250 <- merge(newtab_250, taxALL_250[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_250 = newtab_final_250 %>%
  mutate(Network = "250m")

taxALL_300 <- tax_table(filter_18S_300)
taxALL_300 <- as.data.frame(taxALL_300, check.names=FALSE)
taxALL_300 <- cbind(rownames(taxALL_300), data.frame(taxALL_300, row.names=NULL))
colnames(taxALL_300)[1] <- "Var1"
newtab_300 <- merge(filt_300, taxALL_300[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_300)[1] <- "Var2"
newtab_final_300 <- merge(newtab_300, taxALL_300[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_300 = newtab_final_300 %>%
  mutate(Network = "300m")

taxALL_500 <- tax_table(filter_18S_500)
taxALL_500 <- as.data.frame(taxALL_500, check.names=FALSE)
taxALL_500 <- cbind(rownames(taxALL_500), data.frame(taxALL_500, row.names=NULL))
colnames(taxALL_500)[1] <- "Var1"
newtab_500 <- merge(filt_500, taxALL_500[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_500)[1] <- "Var2"
newtab_final_500 <- merge(newtab_500, taxALL_500[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_500 = newtab_final_500 %>%
  mutate(Network = "500m")

taxALL_600 <- tax_table(filter_18S_600)
taxALL_600 <- as.data.frame(taxALL_600, check.names=FALSE)
taxALL_600 <- cbind(rownames(taxALL_600), data.frame(taxALL_600, row.names=NULL))
colnames(taxALL_600)[1] <- "Var1"
newtab_600 <- merge(filt_600, taxALL_600[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_600)[1] <- "Var2"
newtab_final_600 <- merge(newtab_600, taxALL_600[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_600 = newtab_final_600 %>%
  mutate(Network = "600m")

taxALL_800 <- tax_table(filter_18S_800)
taxALL_800 <- as.data.frame(taxALL_800, check.names=FALSE)
taxALL_800 <- cbind(rownames(taxALL_800), data.frame(taxALL_800, row.names=NULL))
colnames(taxALL_800)[1] <- "Var1"
newtab_800 <- merge(filt_800, taxALL_800[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_800)[1] <- "Var2"
newtab_final_800 <- merge(newtab_800, taxALL_800[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_800 = newtab_final_800 %>%
  mutate(Network = "800m")

taxALL_1000 <- tax_table(filter_18S_1000)
taxALL_1000 <- as.data.frame(taxALL_1000, check.names=FALSE)
taxALL_1000 <- cbind(rownames(taxALL_1000), data.frame(taxALL_1000, row.names=NULL))
colnames(taxALL_1000)[1] <- "Var1"
newtab_1000 <- merge(filt_1000, taxALL_1000[, c("Var1","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_1000)[1] <- "Var2"
newtab_final_1000 <- merge(newtab_1000, taxALL_1000[, c("Var2","Domain","Class", "Order", "Family", "Genus", "Species")], by="Var2")
newtab_final_1000 = newtab_final_1000 %>%
  mutate(Network = "1000m")

# Merge all of the network tables
df_all = rbind(newtab_final_1,newtab_final_40,newtab_final_80,newtab_final_120,newtab_final_160,newtab_final_200,newtab_final_250,newtab_final_300,newtab_final_500,newtab_final_600,newtab_final_800,newtab_final_1000)
df_all$Sign <- ifelse(df_all$Weight > 0, "Pos", "Neg") # Positive sign if edge weight is > 0; if not, it is negative
df_all$Edge = paste(df_all$Class.x, df_all$Class.y,sep="-") # Create edge column for class level

# Subset network edges to Syndiniales and potential host groups
df_syn <- subset(df_all, Edge == "Syndiniales-Arthropoda" | Edge == "Arthropoda-Syndiniales" | Edge == "Syndiniales-Polycystinea" | Edge == "Polycystinea-Syndiniales" | Edge == "Syndiniales-Dinophyceae" | Edge == "Dinophyceae-Syndiniales" | Edge == "Syndiniales-RAD-A" | Edge =="RAD-A-Syndiniales" | Edge == "Syndiniales-RAD-B" | Edge == "RAD-B-Syndiniales" | Edge == "Syndiniales-Acantharea" | Edge == "Acantharea-Syndiniales")
df_syn[df_syn == "Arthropoda-Syndiniales"] <- "Syndiniales-Arthropoda"
df_syn[df_syn == "Polycystinea-Syndiniales"] <- "Syndiniales-Polycystinea"
df_syn[df_syn == "Dinophyceae-Syndiniales"] <- "Syndiniales-Dinophyceae"
df_syn[df_syn == "RAD-A-Syndiniales"] <- "Syndiniales-RAD-A"
df_syn[df_syn == "RAD-B-Syndiniales"] <- "Syndiniales-RAD-B"
df_syn[df_syn == "Acantharea-Syndiniales"] <- "Syndiniales-Acantharea"

# Number of positive/negative edges for each Syndiniales-host pairing at each depth
df_class <- df_syn[c("Class.x", "Class.y")] # Subset to class
sign <- df_syn[c("Sign")] # Subset correlation sign
edge <- df_syn[c("Edge")] # Subset edge
network <- df_syn[c("Network")] # Subset depth network

df_merge = cbind(df_class,edge, sign, network) # Merge pairings, sign, and depth 
df <- as.data.frame(table(df_merge[c(3:5)])) # Subset desired columns
df = df[with(df, order(-Freq)), ] # Estimate the frequency of each type of class level pairing (positive and negative) ordered from largest to smallest

# Plot the number of network edges for each Syndiniales-host pairing with depth - Figure 2C
df$case.weight <- paste0(df$Edge, df$Sign)
df$Network = as.factor(df$Network)
df$Network <- str_remove(df$Network, pattern = "m") # Remove the m after the depth value
p <- ggplot(df, aes(x=fct_rev(Network), y=Freq,color = Edge, lty=Sign,group=case.weight))
p$data$Network <- factor(p$data$Network, levels = c( "1","40","80","120", "160", "200", "250", "300", "500", "600","800","1000"))
p$data$Edge <- factor(p$data$Edge, levels =c("Syndiniales-Arthropoda","Syndiniales-Dinophyceae","Syndiniales-Polycystinea","Syndiniales-Acantharea", "Syndiniales-RAD-A","Syndiniales-RAD-B"))
p + theme_bw() + theme(text = element_text(size=14)) + ylab("# of Network Edges") + xlab("Depth (m)") +theme(legend.position="right")+ 
  scale_color_manual(values=c("#117733", "#CC6677","#DDCC77","#88CCEE","#999933","#882255"))+ 
   geom_point(show.legend = F,size=3)+ geom_line(lwd=1)+ 
  scale_linetype_manual(values = c(Pos = "solid", Neg = "dashed")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, color="black"))+coord_flip() +ylim(0,40) +
  facet_wrap(~Edge)+ theme(legend.key.size = unit(1, 'cm'))
ggsave(filename = "Syn_edges_final.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 10, height = 7, dpi = 300)

# Export the networks and taxonomy files to be analyzed in Cytoscape - we can observe the networks and estimate node degree
write_csv(filt_1, "BATS-network-1m.csv")
write_csv(filt_40, "BATS-network-40m.csv")
write_csv(filt_80, "BATS-network-80m.csv")
write_csv(filt_120, "BATS-network-120m.csv")
write_csv(filt_160, "BATS-network-160m.csv")
write_csv(filt_200, "BATS-network-200m.csv")
write_csv(filt_250, "BATS-network-250m.csv")
write_csv(filt_300, "BATS-network-300m.csv")
write_csv(filt_500, "BATS-network-500m.csv")
write_csv(filt_600, "BATS-network-600m.csv")
write_csv(filt_800, "BATS-network-800m.csv")
write_csv(filt_1000, "BATS-network-1000m.csv")

tax_1 <- as.data.frame(tax_table(filter_18S_surface), check.names=FALSE)
write.table(tax_1,"BATS-network-tax-1m.txt",sep="\t", quote=FALSE)
tax_40 <- as.data.frame(tax_table(filter_18S_40), check.names=FALSE)
write.table(tax_40,"BATS-network-tax-40m.txt",sep="\t", quote=FALSE)
tax_80 <- as.data.frame(tax_table(filter_18S_80), check.names=FALSE)
write.table(tax_80,"BATS-network-tax-80m.txt",sep="\t", quote=FALSE)
tax_120 <- as.data.frame(tax_table(filter_18S_120), check.names=FALSE)
write.table(tax_120,"BATS-network-tax-120m.txt",sep="\t", quote=FALSE)
tax_160 <- as.data.frame(tax_table(filter_18S_160), check.names=FALSE)
write.table(tax_160,"BATS-network-tax-160m.txt",sep="\t", quote=FALSE)
tax_200 <- as.data.frame(tax_table(filter_18S_200), check.names=FALSE)
write.table(tax_200,"BATS-network-tax-200m.txt",sep="\t", quote=FALSE)
tax_250 <- as.data.frame(tax_table(filter_18S_250), check.names=FALSE)
write.table(tax_250,"BATS-network-tax-250m.txt",sep="\t", quote=FALSE)
tax_300 <- as.data.frame(tax_table(filter_18S_300), check.names=FALSE)
write.table(tax_300,"BATS-network-tax-300m.txt",sep="\t", quote=FALSE)
tax_500 <- as.data.frame(tax_table(filter_18S_500), check.names=FALSE)
write.table(tax_500,"BATS-network-tax-500m.txt",sep="\t", quote=FALSE)
tax_600 <- as.data.frame(tax_table(filter_18S_600), check.names=FALSE)
write.table(tax_600,"BATS-network-tax-600m.txt",sep="\t", quote=FALSE)
tax_800 <- as.data.frame(tax_table(filter_18S_800), check.names=FALSE)
write.table(tax_800,"BATS-network-tax-800m.txt",sep="\t", quote=FALSE)
tax_1000 <- as.data.frame(tax_table(filter_18S_1000), check.names=FALSE)
write.table(tax_1000,"BATS-network-tax-1000m.txt",sep="\t", quote=FALSE)

# Upload and plot Syndiniales and host node degree at each depth - Figure 2B and Figure SI 7
major_deg <- read.csv("major_degree.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") 
major_deg$Network = as.factor(major_deg$Network)
df1 <- major_deg %>% # Estimate the mean and SD for each depth and group
  dplyr::group_by(Network,Class) %>%
  dplyr::summarise(
    sd = sd(Degree, na.rm = TRUE),
    len = mean(Degree)
  )
p <- ggplot(df1, aes(x=fct_rev(Network), y=len, fill=Class))
p$data$Class <- factor(p$data$Class, levels =c("Syndiniales","Arthropoda","Dinophyceae","Polycystinea","Acantharea", "RAD-A","RAD-B"))
p + theme_bw() + geom_bar(stat="identity",colour="black")+
  geom_errorbar(aes(ymin = len, ymax = len+sd), width = 0.5) + scale_y_continuous(expand = c(0, 0),limits=c(0,10))+
  theme(text = element_text(size=14)) + ylab("Node degree") + xlab("Depth (m)") +theme(legend.position="right")+ 
  scale_fill_manual(values=c("#332288","#117733", "#CC6677","#DDCC77","#88CCEE","#999933","#882255"))+ 
  facet_wrap(~Class,nrow=2)+coord_flip()
ggsave(filename = "Degree.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 5, dpi = 300)

# Plot total number of network edges connected to Syndiniales at each depth - Figure 2A
edges_tot <- read.csv("edges_bar.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") 
edges_tot$Network = as.factor(edges_tot$Network)
p <- ggplot(data=edges_tot, aes(x=fct_rev(Network), y=Edges, fill=Treat)) 
p$data$Treat <- factor(p$data$Treat, levels = c( "Total","Syn"))
p +geom_bar(stat="identity",colour="black")+scale_y_continuous(expand = c(0, 0),limits=c(0,300))+
  scale_fill_manual(values=c("#8f7bf4","#332288"))+ylab("Network edges") + xlab("Depth (m)")+ theme_bw()+coord_flip()
ggsave(filename = "Total_edges.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 4.5, height = 5, dpi = 300)

# Export the large network edges file (all edges and those connected to Syndiniales) - Table S3
write_csv(df_all, "BATS-network-connections_all.csv") # All network edges
write_csv(df_syn, "BATS-network-connections_syn.csv") # Edges involving Syndiniales and potential hosts

# Re-import file with genus level host dynamics with depth (genera connected to Syndiniales)
genus_inter <- read.csv("BATS_hosts_depth.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") 
df <- as.data.frame(table(genus_inter[c(1:3)])) 
df = df[with(df, order(-Freq)), ] # Estimate frequency of network edges at each depth for genus level 
df[df==0] <- NA # Remove NA values
df<-df[complete.cases(df),]

barplot_hosts = df
barplot_hosts$Genus<- as.character(barplot_hosts$Genus)
barplot_hosts$Network <- str_remove(barplot_hosts$Network, pattern = "m") # Remove m from the depth values
df_arthro <- subset(barplot_hosts, Class =="Arthropoda") # Subset to each respective class level group
df_dino <- subset(barplot_hosts, Class =="Dinophyceae")
df_acanth <- subset(barplot_hosts, Class =="Acantharea")
df_poly <- subset(barplot_hosts, Class =="Polycystinea")
df_rada <- subset(barplot_hosts, Class =="RAD-A")
df_radb <- subset(barplot_hosts, Class =="RAD-B")

# Plot genus level edges with Syndiniales for each group and depth - Figure SI 8
df_arthro$Genus[df_arthro$Freq < 2]= "AAOthers" # Genus level groups that were connected to Syndiniales < 2 times were grouped into "Other" category
p <- ggplot(data=df_arthro, aes(x=fct_rev(Network), y=Freq, fill=Genus))
p$data$Network <- factor(p$data$Network, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000"))
p + geom_bar(aes(), stat="identity", position="stack", width = 0.9)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,50))+ geom_hline(yintercept=0) + theme_bw()+ scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=13, ncol=1)) + coord_flip() +labs(y = "Edge frequency", x = "Depth (m)") 
ggsave(filename = "Syn_genus_arthro.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

df_dino$Genus[df_dino$Freq < 3]= "AAOthers" # Similar plot for Dinophyceae potential hosts
p <- ggplot(data=df_dino, aes(x=fct_rev(Network), y=Freq, fill=Genus))
p$data$Network <- factor(p$data$Network, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000"))
p + geom_bar(aes(), stat="identity", position="stack", width = 0.9)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,50))+ geom_hline(yintercept=0) + theme_bw()+ scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=13, ncol=1)) + coord_flip() +labs(y = "Edge frequency", x = "Depth (m)") 
ggsave(filename = "Syn_genus_dino.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

df_acanth$Genus[df_acanth$Freq < 2]= "AAOthers" # Similar plot for Acantharea potential hosts
p <- ggplot(data=df_acanth, aes(x=fct_rev(Network), y=Freq, fill=Genus))
p$data$Network <- factor(p$data$Network, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000"))
p + geom_bar(aes(), stat="identity", position="stack", width = 0.9)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,50))+ geom_hline(yintercept=0) + theme_bw()+ scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=13, ncol=1)) + coord_flip() +labs(y = "Edge frequency", x = "Depth (m)") 
ggsave(filename = "Syn_genus_acanth.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

df_poly$Genus[df_poly$Freq < 2]= "AAOthers" # Similar plot for Polycystinea potential hosts
p <- ggplot(data=df_poly, aes(x=fct_rev(Network), y=Freq, fill=Genus))
p$data$Network <- factor(p$data$Network, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000"))
p + geom_bar(aes(), stat="identity", position="stack", width = 0.9)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,50))+ geom_hline(yintercept=0) + theme_bw()+ scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=13, ncol=1)) + coord_flip() +labs(y = "Edge frequency", x = "Depth (m)") 
ggsave(filename = "Syn_genus_poly.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

df_rada$Genus[df_rada$Freq < 2]= "AAOthers" # Similar plot for RAD-A potential hosts
p <- ggplot(data=df_rada, aes(x=fct_rev(Network), y=Freq, fill=Genus))
p$data$Network <- factor(p$data$Network, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000"))
p + geom_bar(aes(), stat="identity", position="stack", width = 0.9)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,50))+ geom_hline(yintercept=0) + theme_bw()+ scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=13, ncol=1)) + coord_flip() +labs(y = "Edge frequency", x = "Depth (m)") 
ggsave(filename = "Syn_genus_rada.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

df_radb$Genus[df_radb$Freq < 4]= "AAOthers" # Similar plot for RAD-B potential hosts
p <- ggplot(data=df_radb, aes(x=fct_rev(Network), y=Freq, fill=Genus))
p$data$Network <- factor(p$data$Network, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000"))
p + geom_bar(aes(), stat="identity", position="stack", width = 0.9)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,50))+ geom_hline(yintercept=0) + theme_bw()+ scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1)) + theme(legend.position="right") +  theme(text = element_text(size = 12))+
  guides(fill=guide_legend(nrow=13, ncol=1)) + coord_flip() +labs(y = "Edge frequency", x = "Depth (m)") 
ggsave(filename = "Syn_genus_radb.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

### To prepare for mean POC flux vs. mean group abundance, subset taxonomic data from the photic zone
photic <- x1[x1[["vert_bin"]]== "Photic", ]
mean.18S = photic %>% # Group by class and month (mean and SD)
  dplyr::group_by(Class,Month) %>%
  dplyr::summarize(
    mean = mean(Abundance, na.rm = FALSE), 
    sd = sd(Abundance, na.rm = FALSE))

# Subset to Syndiniales, Arthropoda, and Dinophyceae
abund_18S <- subset(mean.18S, grepl("Syndiniales|Dinophyceae|Arthropoda", mean.18S$Class)) 
write_csv(abund_18S, "BATS-vert-groups.csv") # Save file for mean group relative abundance in the photic zone

# Import mean relative abundances and POC flux data (150 m) and plot Spearman correlation - Figure 3
syn_export <- read.csv("Syn_export.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") 
p <- ggscatter(syn_export, x = "C_avg", y = "Mean", cor.method = "spearman", cor.coef =T, size = 3, add.params = list(color = "black", fill = "darkgray"), add = "reg.line", conf.int = TRUE, ylab = "Integrated relative abundance <140 m (%)", xlab = "Carbon flux 150 m (mgC/m^2/day)")+
  geom_errorbar(data = syn_export, aes(ymin = Mean-sd, ymax = Mean+sd),width = 2)+ # Error bars based on average relative abundance and SD
  geom_errorbar(data = syn_export, aes(xmin = C_avg-C_sd, xmax = C_avg+C_sd),width = 2)+ # Error bars based on average POC flux and SD
  geom_point(aes(fill=Class),size = 5, shape=21,colour = "black")  + theme_bw() + theme(text = element_text(size = 14)) +facet_grid(~factor(Class, levels=c('Syndiniales', 'Dinophyceae', 'Arthropoda')))+
  scale_fill_manual(values=c("#332288","#117733","#CC6677"))
p$data$Class <- factor(p$data$Class, levels =c("Syndiniales","Arthropoda","Dinophyceae"))     
p
ggsave(filename = "EXPORT_all.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 12, height = 4, dpi = 300)

# Merge phyloseq objects from the network to run an UpSet plot - overlap of Syndiniales presence-absence data with depth - Figure 4A
ps_all = merge_phyloseq(filter_18S_surface,filter_18S_40,filter_18S_80,filter_18S_120,filter_18S_160,filter_18S_200,filter_18S_250,filter_18S_300,filter_18S_500,filter_18S_600,filter_18S_800,filter_18S_1000)
ps_Syn_all =subset_taxa(ps_all, Class=="Syndiniales",Prune = T) # Subset to Syndiniales
dataset <- phyloseq2meco(ps_Syn_all)
dataset1 <- dataset$merge_samples(use_group = "Nominal_Depth")
otu <- as.data.frame(dataset1$otu_table)
otu2 = otu[,c(2,12,10,9,7,6,5,4,3,11,8,1)]
depth = colnames(otu2)
pdf("Syn_upset.pdf", width = 15, height = 8)
upset(otu,depth, width_ratio=0.1,min_size=1, sort_sets = FALSE,name = "Syndiniales ASVs",base_annotations=list('Size'=(intersection_size(counts=FALSE))),queries=list(
  upset_query(set='1', fill="#DDAA33",color='black'),
  upset_query(set='40', fill="#DDAA33",color='black'),
  upset_query(set='80', fill="#DDAA33",color='black'),
  upset_query(set='120', fill="#DDAA33",color='black'),
  upset_query(set='160', fill="#004488",color='black'),
  upset_query(set='200', fill="#004488",color='black'),
  upset_query(set='250', fill="#004488",color='black'),
  upset_query(set='300', fill="#004488",color='black'),
  upset_query(set='500', fill="#004488",color='black'),
  upset_query(set='600', fill="#004488",color='black'),
  upset_query(set='800', fill="#004488",color='black'),
  upset_query(set='1000', fill="#004488",color='black')
),  matrix=intersection_matrix(
  geom=geom_point(
    shape='circle filled',
    size=5,
    stroke=0.5
  ),
))

dev.off()

# Read counts of two example Syndiniales ASVs over time and depth at BATS - Figure 4B
ps_melt2 <- psmelt(ps_Syn_all) # Melt the data
ps_melt3 <- subset(ps_melt2, OTU == "ASV9865" | OTU == "ASV2889") # Subset to two specific ASVs - chosen from the UpSet plot
p <- ggplot(ps_melt3, aes(x=fct_rev(Nominal_Depth), y = factor(yyyymmdd),fill=vert_bin))
p$data$OTU <- factor(p$data$OTU, levels = c("ASV2889","ASV9865")) # Order the ASVs on plot
p$data$Nominal_Depth <- factor(p$data$Nominal_Depth, levels = c( "1", "40", "80","120","160","200","250","300","500","600","800","1000")) # Order depths
p + geom_point(aes(size=Abundance, fill= vert_bin), alpha = 1, shape = 21) + theme_bw() + scale_fill_manual(values=c("#004488","#DDAA33"))+
  theme(axis.text.x=element_text(angle=45,vjust =1, hjust=1, size=12)) +theme(axis.text.y=element_text(size =12))+
  scale_size_continuous(limits = c(1, 2500), range=c(1,6), breaks = c(1, 500,1000, 1500,2000,2500)) +
  theme(axis.title.x = element_blank())+labs(x="Depth (m)")+
  guides(fill=guide_legend(ncol=2, override.aes = list(size=3))) + coord_flip() +
  facet_wrap(~OTU,nrow=2) 
ggsave(filename = "Syn_ASVs_export.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 7, height = 5, dpi = 300)

# Estimate the percent of Syndiniales, and other groups, that are present in the photic and aphotic zones
dataset <- phyloseq2meco(ps_rare)
dataset$tax_table %<>% .[grepl("c__Syndiniales", .$Class), ] # Subset to Syndiniales  
dataset$tax_table %<>% .[grepl("c__Arthropoda", .$Class), ] # Subset to Arthropoda  
dataset$tax_table %<>% .[grepl("c__Dinophyceae", .$Class), ] # Subset to Dinophyceae
dataset1 <- dataset$merge_samples(use_group = "vert_bin") # Partition based on vertical zones
t1 <- trans_venn$new(dataset1, ratio = "numratio") # Apply venn diagram - move through different groups
t1$plot_venn(color_circle = c("#88CCEE","#CC6677"),fill_color = TRUE,alpha=0.5) # Plot overlap

# Upload overlap data for each group and plot - Figure 4C
shared <- read.csv("shared_ASVs.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") 
p <- ggplot(shared, aes(x=Group, y=as.numeric(Shared),fill=Group))
p$data$Group <- factor(p$data$Group, levels = c("All euks","Syndiniales","Dinophyceae","Arthropoda"))
p + geom_col(width = 0.8, colour = "black", size = 0.9)  + theme_bw()+scale_fill_manual(values=c("#757575","#332288","#CC6677","#117733","#DDCC77","#88CCEE","#999933","#882255"))+ geom_hline(yintercept=20, linetype=2, colour="darkred")+
  scale_y_continuous(name="% of ASVs present > and < 200 m", limits=c(0, 30), expand = c(0, 0)) + theme(axis.title.x=element_blank()) + 
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black"))+ theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 12))
ggsave(filename = "shared_ASVs.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 300)

# UpSet plot to look at seasonal effects - Figure SI 9
dataset <- phyloseq2meco(ps_Syn_all)
dataset1 <- dataset$merge_samples(use_group = "Season")
otu <- as.data.frame(dataset1$otu_table)
otu2 = otu[,c(1:4)]
depth = c("Fall","Summer","Spring","Winter") # Set order of seasons
upset(otu,depth, width_ratio=0.1,min_size=1, sort_sets = FALSE,name = "Syndiniales ASVs",base_annotations=list('Size'=(intersection_size(counts=TRUE))),queries=list(
  upset_query(set='Winter', fill=ap[5],color=ap[5]),
  upset_query(set='Spring', fill=ap[4],color=ap[4]),
  upset_query(set='Summer', fill=ap[12],color=ap[12]),
  upset_query(set='Fall', fill=ap[8],color=ap[8])
),  matrix=intersection_matrix(
  geom=geom_point(
    shape='circle filled',
    size=3.5,
    stroke=0.45
  ),
))

# Run PCoA at select depths at BATS and partition samples based on season and year - Figure SI 3
depth = subset_samples(ps_Syn, Nominal_Depth=="1") # Repeat for 120, 600 and 1000 m
ordu = ordinate(depth, "PCoA", "bray")
p = plot_ordination(depth, ordu, color="Season",shape="Year")
p$data$Year <- as.factor(p$data$Year)
p$data$Season <- as.factor(p$data$Season)
p + theme_bw() + 
  scale_colour_manual(values=c("#DB5339", "#A5CFCC", "#F7CDA4", "#0E899F")) + scale_shape_manual(values = c(15:18))+
  geom_point(size=6) + theme(text = element_text(size=12)) 
ggsave(filename = "Syn_PCoA_1m.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 7, height = 5, dpi = 300)

# PERMANOVA at different depths - temporal component
metadata <- as(sample_data(depth), "data.frame")
metadata$Season=as.factor(metadata$Season)
adonis2(phyloseq::distance(depth, method="bray") ~ Season, data = metadata,perm=9999)

# Re-run taxonomic bar plots of the whole community and Syndiniales only filtered to the top 150 ASVs from each depth used in network analysis 
# Plot the full community used in the networks - Figure SI 6
dataset <- phyloseq2meco(ps_all) # Convert to microeco object
dataset1 <- dataset$merge_samples(use_group = "Nominal_Depth") # Group the data by depth
t1 <- trans_abund$new(dataset = dataset1, taxrank = "Class", ntaxa = 10) # Create an abundance object with the top 10 class level 18S groups
t1$data_taxanames <- c("Syndiniales", "Dinophyceae", "Polycystinea","Arthropoda", 
                       "Acantharea", "RAD-B","Cnidaria","RAD-A","Sagenista","Prymnesiophyceae")
g1 <- t1$plot_bar(bar_type = "full", use_alluvium = FALSE, clustering = FALSE,barwidth = NULL, xtext_size = 12, color_values = c("#332288","#CC6677","#DDCC77","#117733","#88CCEE", "#882255", "#44AA99", "#999933","#AA4499" , "#EB7357"),others_color = "#757575",
                  order_x = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))
g1 + coord_flip() + geom_col(colour = "black") + labs(x="Depth (m)",y="Relative abundance (%)")+theme(text = element_text(size = 12))
ggsave(filename = "syn_stacked_class_network.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 5, dpi = 300)

# Plot Syndiniales used in the networks - Figure SI 6
dataset <- phyloseq2meco(ps_Syn_all) # Convert to microeco object
dataset1 <- dataset$merge_samples(use_group = "Nominal_Depth") # Group the data by depth
t1 <- trans_abund$new(dataset = dataset1, taxrank = "Family", ntaxa = 10) 
t1$data_taxanames <- c("Dino-Group-II-Clade-7", "Dino-Group-II-Clade-6", "Dino-Group-II-Clade-10-and-11","Dino-Group-I-Clade-1", 
                       "Dino-Group-I-Clade-5", "Dino-Group-I-Clade-2","Dino-Group-II_X","Dino-Group-I-Clade-4","Dino-Group-I_X","Dino-Group-II-Clade-22")
g1 <- t1$plot_bar(bar_type = "full", use_alluvium =FALSE, clustering = FALSE,barwidth = NULL, xtext_size = 12, color_values = c("#4477AA","#EE6677","#228833", "#CCBB44", "#66CCEE", "#AA3377", "#EE7733","#CC3311","#009988","#EE3377"),others_color = "#757575",
                  order_x = c("1000", "800", "600", "500", "300", "250", "200","160","120","80","40","1"))
g1 + coord_flip()+geom_col(colour = "black")+ labs(x="Depth (m)",y="Relative abundance (%)")+theme(text = element_text(size = 12))
ggsave(filename = "syn_stacked_family_network.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5.5, height = 5, dpi = 300)
