#############################################################################
#      Microbiome analysis for Daphnia Toxic vs. Non Toxic Exposure         #
#    Last Updated 8/11/2024 by Amir Gabidulin (amir.gabidulin@wsu.edu)      #
# See associated manuscript for more details about the analysis and results #
#############################################################################
#Essentials Packages:
library(tidyverse)
library(dplyr)
library(vegan)
library(phyloseq)
library(speedyseq)
library(ape)
library(readxl)
library(lme4)
library(lmerTest)
library(DESeq2)
library(ggpubr)

#Miscellaneous Packages (Not required, but load if some errors occur):
library(phylosmith)
library(ggpubr)
library(visreg)
library(car)
library(ResourceSelection)
library(mlmRev)
library(report)
library(reshape2)
library(hrbrthemes)
library(viridis)
library(extrafont)
library(extrafontdb)
library(ggplot2)
library(emmeans)
library(forecast)
library(colorspace)
library(corrr)
library(cowplot)
library(ggdark)
library(ggforce)
library(ggrepel)
library(ggridges)
library(ggsci)
library(ggtext)
library(ggthemes)
library(grid)
library(gridExtra)
library(patchwork)
library(scico)
library(showtext)
library(shiny)
library(psych)
library(rstatix)
library(correlation)
library(palmerpenguins)
library(wesanderson)
library(RColorBrewer)
library(ggthemes)
library(paletteer)
library(cartography)
library(qiime2R)
library(ggalt)
library(cxhull)
library(ggforce)
library(tibble)
library(tidyr)
library(stringr)
library(remotes)
library(dada2)
library(MicEco)
library(DBI)
library(mia)
library(metagenomeSeq)
library(ANCOMBC)
library(ALDEx2)
library(microbiomeMarker)
library(proxy)


#### ____ Importing Genus, OTU Table, Metadata, and Phylogenetic Tree ____ ####:
#Phylogeny 2: 

import_taxa <- read.table('your/directory/taxonomy.tsv', header = TRUE, sep = '\t')
head(import_taxa)

# First we have to provide names for the new columns
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxonomy <- import_taxa %>%
  mutate_at('Taxon',str_replace_all, "[a-z]__","") %>%
  separate(Taxon, sep = ';', into=ranks,remove = TRUE) %>%
  column_to_rownames(var = "Feature.ID") %>%
  as.matrix()
head(taxonomy)

TAX = tax_table(taxonomy)

#OTU Table import 2 (MAIN): 
import_table <- read.table('your/directory/phyloseq-data/feature-table.txt', header = TRUE, sep = '\t', row.names = 1, comment.char = "")
otumat <- as.matrix(import_table)
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)

#Load Metadata:
metadata <- read.table("your/directory/DaphniaMetaData123.txt", header = T, sep = '\t')

#Read the metadata file into phyloseq:
metadata[] <- lapply(metadata, as.factor)
head(metadata)
row.names(metadata) <- paste0('X', gsub('-','.', metadata$SampleID))
tail(metadata)
META <- sample_data(metadata)

# Load Rooted Phylogenetic Tree:
otu_tree <- read.tree(file = "your/directory/phyloseq-data/tree.nwk")
plot(otu_tree)

# Create Phyloseq Object:
physeq <- phyloseq(OTU,TAX,META,otu_tree)
physeq <- prune_taxa(taxa_sums(physeq) > 5, physeq)
physeq

#Combine Replicates across metadata column "Seperation":
physeq_combined <- merge_samples2(physeq, "Seperation", fun_otu = sum)
sample_data(physeq_combined)
physeq_combined

#### ____ Rarefaction assessment ____ ####:
#Generate rarefaction quantiles:
qunatiles <- seq(0,1,0.05)
percentiles <- quantile(sample_sums(physeq_combined), qunatiles)
print(percentiles)

#Species accumulation curve:
otu_table <- otu_table(physeq_combined)
otu_data <- as.data.frame(otu_table)
total_depth <- rowSums(otu_data)
spec_accum <- specaccum(otu_data)
# Plot the rarefaction curve
plot(spec_accum, col = "blue", xlab = "Number of Samples", ylab = "Observed OTUs",
     main = "Rarefaction Curve")

#Rarefication Curve:
tab <- otu_table(physeq_combined)
class(tab) <- "matrix"
tab <- t(tab)
rare <- rarecurve(tab, step = 100, ylab="ASV", label=F)

#Find what data is below a certain threshold:
sample_sums_values <- sample_sums(physeq_combined)
threshold <- 20000
below_threshold_ids <- names(sample_sums_values[sample_sums_values < threshold])
print(below_threshold_ids)

#Rarefy phyloseq object:
physeq.rarefied <- rarefy_even_depth(physeq_combined, rngseed=1, sample.size= 20000)
plot_bar(physeq.rarefied)
physeq.rarefied
sample_data(physeq.rarefied)


#### ____ Plotting Relative Abundances ____ ####:
ps = physeq.rarefied
#ps = prune_taxa(taxa_sums(ps) > 1000, ps)

total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_ps = transform_sample_counts(ps, standf)

normalized_ps_phylum1 = tax_glom(normalized_ps, "Class")
normalized_ps_phylum1
normalized_ps_phylum1_relabun = transform_sample_counts(normalized_ps_phylum1, function(OTU) OTU/sum(OTU) * 100)

library(reshape2)
cud_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#E69F00", 
  "#56B4E9", "#009E73", "#F0E442", "#0072B2"
)

#Abundances with facet_wrap(~Treatment)
psmelt(normalized_ps_phylum1_relabun) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cud_palette) +  # Using the CUD palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Relative Abundance", title = "Class Relative Abundance")

# Filter out Class level objects that are under 1% of the whole representation
threshold <- 1  # set your threshold percentage here

filtered_ps <- prune_taxa(taxa_sums(normalized_ps_phylum1_relabun) > threshold, normalized_ps_phylum1_relabun)

CData <- subset_samples(filtered_ps, Treatment=="C")
CData

GData <- subset_samples(filtered_ps, Treatment=="G")
GData


# Plotting Relative Abundances C Treatment
graph1 <- psmelt(CData) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cud_palette) + theme_bw() +# Using the CUD palette
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(y = "Relative Abundance (%)") + facet_wrap(~Treatment)

# Plotting Relative Abundances G Treatment
graph2 <- psmelt(GData) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cud_palette) + theme_bw() +  # Using the CUD palette
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank()) +
  labs(y = "Relative Abundance (%)") + facet_wrap(~Treatment) 

# Combine Relative Abundance Graphs by Treatment
abgraph <- ggarrange(graph1, graph2, common.legend = TRUE, legend = "right")
abgraph



####################################################
#### ____ Alpha and Beta Diversity Metrics ____ ####:
####################################################
####Alpa Diversity Metics: 
alpha_shannon <- estimate_richness(physeq.rarefied, measures = "Shannon")
alpha_analysis <- lm(alpha_shannon$Shannon ~ Clone, data = metaData123)
alpha_analysis <- lm(alpha_shannon$Shannon ~ Treatment, data = metaData123)
summary(alpha_analysis)

####Bray_Curtis:
distance_matrix_bray <- phyloseq::distance(physeq.rarefied, method = "bray")
#Ordination Alignment:
GP.ord <- ordinate(physeq.rarefied, "PCoA", "bray")
p1 = plot_ordination(physeq.rarefied, GP.ord, type = "taxa", color = "Genus", title = "taxa")
p1 + facet_wrap(~Phylum, 3)
print(p1)
#Metadata Alignment:
GP.ord <- ordinate(physeq.rarefied, "PCoA", "bray")
p1 = plot_ordination(physeq.rarefied, GP.ord, label = "Clone", type = "samples", color = "Treatment", title = "Treatment Exposure PCoA(Bray)")
p1 + facet_wrap(~Clone)
print(p1 + theme_bw())


####Weighted_Unifrac:
#Generate Distance Matrix; 
distance_matrix_unifrac_weighted <- phyloseq::distance(physeq.rarefied, method = "unifrac", weighted = TRUE)
distance_matrix_bray <- phyloseq::distance(physeq.rarefied, method = "bray")
#Oridnation:
unifrac.ord <- ordinate(physeq.rarefied, "PCoA", "unifrac", weighted = TRUE)
p3 = plot_ordination(physeq.rarefied, unifrac.ord, type = "taxa", color = "Genus", title = "taxa")
p3 + facet_wrap(~phylum, 3)
print(p3)
#Metadata:
unifrac.ord <- ordinate(physeq.rarefied, "PCoA", "unifrac", weighted = TRUE)
p3 = plot_ordination(physeq.rarefied, unifrac.ord, label = "Clone" , type = "samples", color = "Treatment", title = "Treatment Exposure PCoA (Weighted-Unif)")
p3 + facet_wrap(~Clone)
print(p3  + theme_bw())



####Permanovas: 
metaData123 <- sample_data(physeq.rarefied)
metaData123 <- as.matrix.data.frame(metaData123)
metaData123 <- as.data.frame(metaData123)
metaData123$NumberofGuts <- as.numeric(metaData123$NumberofGuts)
metaData123$Clone <- as.numeric(metaData123$Clone)

weighted_unifrac <- adonis2(distance_matrix_unifrac_weighted ~ Treatment + Clone + (1|NumberofGuts), data = metaData123, permutations = 1000)
print(weighted_unifrac)

bray_pr <- adonis2(distance_matrix_bray  ~ Treatment + Clone + (1|NumberofGuts), data = metaData123, permutations = 1000)
print(bray_pr)


unifrac.ord <- ordinate(physeq.rarefied, "PCoA", "unifrac", weighted = TRUE)
GP.ord <- ordinate(physeq.rarefied, "PCoA", "bray")


####Multivariate parallelism plasticity analysis:
unifrac.ord$vectors

#Betadispers
mod <- betadisper(distance_matrix_unifrac_weighted, metaData123$Treatment)
permutest(mod)

mod <- betadisper(distance_matrix_bray, metaData123$Treatment)
permutest(mod)

####Harrar Script Input Preparation:
View(sample_data(physeq.rarefied))

#Generate PCA using distance matrix (weighted-unifrac) (Haerar Script):
wu_pcoa <- as.data.frame(unifrac.ord$vectors)
wu_pcoa_metadata <- sample_data(physeq.rarefied)
wu_pcoa_merged <- merge(wu_pcoa, wu_pcoa_metadata, by = "row.names", all = TRUE)
wu_pcoa_merged <- wu_pcoa_merged %>%
  group_by(Clone) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(group = paste0(gsub("_", "", Treatment), "_", Clone)) %>%
  rename_with(~ gsub("Axis\\.", "PC", .), starts_with("Axis."))
View(wu_pcoa_merged)
write.csv(wu_pcoa_merged, "WU_HarrarSciprt-Daphnia.csv")
#Modify the Output according to Harrar example data


#Generate PCA using distance matrix (Bray) (Haerar Script):
GP_pcoa <- as.data.frame(GP.ord$vectors)
GP_pcoa_metadata <- sample_data(physeq.rarefied)
GP_pcoa_merged <- merge(GP_pcoa, GP_pcoa_metadata, by = "row.names", all = TRUE)
GP_pcoa_merged <- GP_pcoa_merged %>%
  group_by(Clone) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(group = paste0(gsub("_", "", Treatment), "_", Clone)) %>%
  rename_with(~ gsub("Axis\\.", "PC", .), starts_with("Axis."))
View(GP_pcoa_merged)
write.csv(GP_pcoa_merged, "GP_HarrarSciprt-Daphnia.csv")
#Modify the Output according to Harrar example data

##########################################################################################################################################
#Harrar Script (Modified by Amir Gabidulin - Original Script did not perform correctly on large datasets (greater than 20 replicates))   #:
#                 Read in the function first (starts with #theta function & ends with multivariate.vector function)                      #
##########################################################################################################################################
library("dplyr")
library("plyr")
library("tibble")
library("data.table")

dg <- read.csv("your/directory/WU_HarrarSciprt-Daphnia.csv")
dg1 <- read.csv("your/directory/GP_HarrarSciprt-Daphnia.csv")


#Weighted Unifrac Parallelism Testing:
plasticityUW <- multivariate.vector(dg)
WUAngle <- plasticityUW$theta
row_sums <- rowSums(WUAngle)
WU_Angle_result <- row_sums / (ncol(WUAngle) - 1)
t.test(WU_Angle_result, mu = 90)

#Bray Curtis Parallelism Testing:
plasticityGP <- multivariate.vector(dg1)
GPAngle <- plasticityGP$theta
row_sums <- rowSums(GPAngle)
GP_Angle_result <- row_sums / (ncol(GPAngle) - 1)
t.test(GP_Angle_result, mu = 90)

#### ____↓↓↓↓ Read this in first ↓↓↓↓____ ####:
#theta function
f.theta <- function(table.of.diff.in.means, unique.id, select.col){
  angle.matrix.radians <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
      print(i)
    for(j in 1:length(unique.id)){
      angle.matrix.radians[i,j] <- round(acos(cor(x = t(table.of.diff.in.means[i,select.col]), t(table.of.diff.in.means[j,select.col]), use = "na.or.complete")), 3) #
      #angle.practice <- cor(x = t(table.of.diff.in.means[i,select.col]), t(table.of.diff.in.means[j,select.col]), use = "na.or.complete")
    }
  }
  
  rownames(angle.matrix.radians) <- unique.id
  colnames(angle.matrix.radians) <- unique.id
  angle.matrix.degrees <- round(angle.matrix.radians*(180/pi), 3)
  angle.output <- list(angle.matrix.radians, angle.matrix.degrees)
  names(angle.output) <- c("theta.radians", "theta.degrees")
  return(angle.output)
}

#meanL function
f.meanL <- function(table.of.diff.in.means, unique.id, select.col){
  get.vectlength <- function(vect){
    return(sqrt(sum(vect^2, na.rm = TRUE)))
  }
  table.of.diff.in.means.t <- table.of.diff.in.means[,select.col]
  length.of.vector.by.group <- apply(table.of.diff.in.means.t, 1 ,get.vectlength) #1 signifies by row. apply the function to every row in the input dataframe
  length.diff.matrix <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
    for(j in 1:length(unique.id)){
      length.diff.matrix[i,j] <-mean(round(c(length.of.vector.by.group[i], length.of.vector.by.group[j]), 3))
    }
  }
  rownames(length.diff.matrix) <- unique.id
  colnames(length.diff.matrix) <- unique.id
  
  length.of.vector.by.group <- as.data.frame(length.of.vector.by.group)
  length.of.vector.by.group$wshd <- rep(unique.id, length.out = nrow(length.of.vector.by.group))
  #length.of.vector.by.group$wshd <- unique.id
  length.output <- list(length.diff.matrix, length.of.vector.by.group)
  return(length.output)
}

diff <- function(x) {x-lag(x)}

#Vector analysis function
multivariate.vector <- function(infile){
  
  infile2 <- infile[,-c(1:3)]
  
  infile_mean <- aggregate(infile2[,2:length(colnames(infile2))], list (group = infile2$group), mean)
  
  infile_mean2 <- infile_mean %>% add_column(col1 = NA, .after = "group")
  infile_mean2 <- infile_mean2 %>% add_column(col2 = NA, .after = "group")
  
  colnames(infile_mean2)[2] <- colnames(infile)[1]
  colnames(infile_mean2)[3] <- colnames(infile)[2]
  
  infile_mean2[,2] <- sapply(strsplit(infile_mean2$group,"_"), `[`, 1)
  infile_mean2[,3] <- sapply(strsplit(infile_mean2$group,"_"), `[`, 2)
  
  infile_mean2[,1] <- as.factor(infile_mean2[,1])
  infile_mean2[,2] <- as.factor(infile_mean2[,2])
  infile_mean2[,3] <- as.factor(infile_mean2[,3])
  
  infile_mean3 <- infile_mean2[,-c(1,3)]
  
  pcoas <- colnames(infile_mean3[2:length(colnames(infile_mean3))])
  
  differences <- infile_mean3 %>% group_by(infile_mean3[1]) %>% mutate_at(pcoas, function(x) (x- shift(x)))
  
  differences <- as.data.frame(differences)
  differences2 <- differences[,-1]
  i <- (rowSums(differences2,na.rm=T) !=0)
  differences3 <- differences2[i,]
  
  unique.id <- c(unique(infile_mean2[,2]))
  
  #estimate theta
  table.of.diff.in.means <- differences3
  
  select.col <- c(1:length(colnames(differences3)))
  x.theta <- f.theta(table.of.diff.in.means, unique.id, select.col)
  
  cat("Angles between vectors\n")
  print(x.theta$theta.degrees)
  
  cat("\n")
  
  #estimate mean length
  x.mean_length <- f.meanL(table.of.diff.in.means, unique.id, select.col)
  cat("Mean vector lengths\n")
  print(x.mean_length[[1]])
  
  # Create a list to hold the dataframes:
  results_list <- list(theta = x.theta$theta.degrees, length = x.mean_length[[1]])
  return(results_list)
}
#### ____↑↑↑↑ Read this in first ↑↑↑↑____ ####:


########################################
#DESeq2 differential abundance analysis#:
########################################
alpha <- 0.01

diagdds <- phyloseq_to_deseq2(physeq.rarefied, ~ Treatment)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")


res <- results(diagdds, cooksCutoff = FALSE)
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as (tax_table(physeq.rarefied)[rownames(sigtab), ], "matrix"))
#sigtab <- cbind(as(sigtab, "data.frame"), as (sample_data(physeq.rarefied)[rownames(sigtab), ], "matrix"))
View(sigtab)
sigtab$pvalue <- p.adjust(sigtab$pvalue, method = "BH")

sigtab$mean_norm_count <- (sigtab$baseMean + 1) * sigtab$log2FoldChange

diffabb <- ggplot(sigtab, aes(x = mean_norm_count, y = log2FoldChange, color = Family)) +
  geom_point(alpha = 1, size = 3) + scale_fill_manual(values = cud_palette) +
  labs(x = "Mean Normalized Count", y = "Log2 Fold Change") +
  theme_classic()# + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

diffabb

# Set the dimensions for the output image
png_width <- 3000
png_height <- 4500
800
1200

# Adjust the size of the plots in ggarrange using the widths parameter
ggarrange(graph1, graph2, diffabb, labels = c("A","","B"),
          heights = c(1),
          align = "hv", nrow = 2, ncol = 2, scale = 0.5)  # Adjust the width of the second plot here

ggsave("arranged_plots.png", width = png_width, height = png_height)


ggarrange(
  ggarrange(graph1, graph2, 
            nrow = 1, ncol = 2, common.legend = TRUE, legend = 'right'),  # Arrange graph1 and graph2 with a common legend
  diffabb + theme(legend.position = "right"),  # Arrange diffabb with its own legend
  nrow = 2, ncol = 1,
  heights = c(1, 0.5) ,  # Set relative heights for graph1, graph2, and diffabb
  align = "hv"
)


graph1 <- graph1 +ggtitle("A")
graph2 <- graph2 +ggtitle("")
diffabb <- diffabb + ggtitle("B")

top_row <- ggarrange(graph1,graph2, common.legend = TRUE, legend = 'right')
arranged_plot <- ggarrange(top_row, diffabb, ncol = 1, nrow=2)

ggsave("Figure3.png", arranged_plot, width = 10, height = 15, dpi = 300)

#################################
#R0 and Distance Matrix Analysis#:
#################################
R01 <- read.csv("DaphniaMicrobiomeNetReproductiveRate_240228.csv")
R01$treatment <- as.factor(R01$treatment)
hist(R01$R0.Delta)
r0lm <- lm(R0.Delta ~ clone + treatment, data = R01)
summary(r0lm)

diag_values <- diag(plasticityGP$length)
diag_values


daphniaplast <- readxl::read_xlsx("your/directory/DaphniaRNaught.xlsx")
head(daphniaplast)
daphniaplast$clone <- as.factor(daphniaplast$clone)

#UW
lm_model <- lm(`Mean-Vector-Length-UW` ~ `R0-Delta`, data = daphniaplast)
summary(lm_model)
r_squared <- summary(lm_model)$r.squared
r_squared

p <- ggplot(daphniaplast, aes(x = `R0-Delta`, y = `Mean-Vector-Length-UW`)) +
  geom_point() +  # Scatter plot of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model fit line
  labs(
    title = "Relationship Between R0-Delta and Mean-Vector-Length-UW",
    x = "R0-Delta",
    y = "Mean-Vector-Length-UW"
  )
p1 <- p + annotate(
  "text",
  x = max(daphniaplast$`R0-Delta`, na.rm = TRUE),
  y = min(daphniaplast$`Mean-Vector-Length-UW`, na.rm = TRUE),
  label = paste("R-squared =", round(r_squared, 3)),
  hjust = 1,
  vjust = 0
)
p1

#BC
lm_modelBC <- lm(`Mean-Vector-Length-BC` ~ `R0-Delta`, data = daphniaplast)
summary(lm_modelBC)
r_squared <- summary(lm_modelBC)$r.squared
r_squared
p <- ggplot(daphniaplast, aes(x = `R0-Delta`, y = `Mean-Vector-Length-BC`)) +
  geom_point() +  # Scatter plot of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model fit line
  labs(
    title = "Relationship Between R0-Delta and Mean-Vector-Length-BC",
    x = "R0-Delta",
    y = "Mean-Vector-Length-BC"
  )
p2 <- p + annotate(
  "text",
  x = max(daphniaplast$`R0-Delta`, na.rm = TRUE),
  y = min(daphniaplast$`Mean-Vector-Length-BC`, na.rm = TRUE),
  label = paste("R-squared =", round(r_squared, 3)),
  hjust = 1,
  vjust = 0
)
p2
############################
#Neonate Product instead R0#:
############################
#UW
lm_model <- lm(`Mean-Vector-Length-UW` ~ `Mean_NeonateDiff`, data = daphniaplast)
summary(lm_model)
r_squared <- summary(lm_model)$r.squared
r_squared

p <- ggplot(daphniaplast, aes(x = `Mean_NeonateDiff`, y = `Mean-Vector-Length-UW`)) +
  geom_point() +  # Scatter plot of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model fit line
  labs(
    title = "Relationship Between Mean_NeonateDiff and Mean-Vector-Length-UW",
    x = "Mean_NeonateDiff",
    y = "Mean-Vector-Length-UW"
  )
p + annotate(
  "text",
  x = max(daphniaplast$`Mean_NeonateDiff`, na.rm = TRUE),
  y = min(daphniaplast$`Mean-Vector-Length-UW`, na.rm = TRUE),
  label = paste("R-squared =", round(r_squared, 3)),
  hjust = 1,
  vjust = 0
)

#BC
lm_modelBC <- lm(`Mean-Vector-Length-BC` ~ `Mean_NeonateDiff`, data = daphniaplast)
summary(lm_modelBC)
r_squared <- summary(lm_modelBC)$r.squared
r_squared
p <- ggplot(daphniaplast, aes(x = `Mean_NeonateDiff`, y = `Mean-Vector-Length-BC`)) +
  geom_point() +  # Scatter plot of the points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model fit line
  labs(
    title = "Relationship Between Mean_NeonateDiff and Mean-Vector-Length-BC",
    x = "Mean_NeonateDiff",
    y = "Mean-Vector-Length-BC"
  )

p + annotate(
  "text",
  x = max(daphniaplast$`Mean_NeonateDiff`, na.rm = TRUE),
  y = min(daphniaplast$`Mean-Vector-Length-BC`, na.rm = TRUE),
  label = paste("R-squared =", round(r_squared, 3)),
  hjust = 1,
  vjust = 0
)



######################### PhyloFactor #########################
library(phylofactor)
library(phytools)


#Importing Genus, OTU Table, Metadata, and Phylogenetic Tree:
#Phylogeny 2: 

OTU_factorization <- otu_table(physeq.rarefied)
pf_tree <- phy_tree(physeq.rarefied)
mapfile <- as.data.frame(sample_data(physeq.rarefied))
mapfile <- mapfile$Treatment
#pf_tax <- tax_table(physeq.rarefied)
#TAX

pf_tax <- read.table('your/directory/phyloseq-data/taxonomy.tsv', header = TRUE, sep = '\t')
pf_tax <- subset(pf_tax, select = -Confidence)

PF <- PhyloFactor(OTU_factorization, pf_tree, mapfile_clean, nfactors=12, frmla = Data~Treatment)
smry <- pf.summary(PF, factor=1, pf_tax)
pf.tidy(smry)
PF
names(smry)

par(mfrow=c(1,2))
plot(mapfile,smry$ilr,main='Isometric log-ratio',ylab='ILR balance')
plot(mapfile,smry$MeanRatio,main='Ratio of Geometric Means',ylab='Group1/Group2')


smry$TaxaSplit %>%
  lapply(.,FUN=function(x) unique(x$TaxaIDs))

Group1.otus <- smry$group1$IDs$otuIDs %>% as.list %>% sapply(.,toString)
Group1.otus

filtered_pf_tax <- pf_tax %>%
  filter(Feature.ID %in% Group1.otus)
View(filtered_pf_tax)


findElbow <- function(y, plot = FALSE, returnIndex = TRUE) {
  
  # The following helper functions were found at
  # paulbourke.net/geometry/pointlineplane/pointline.r
  # via the SO reference below.  The short segment check
  # was modified to permit vectorization.
  
  ##========================================================
  ##
  ##  Credits:
  ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
  ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
  ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
  ##  With grateful thanks for answering our needs!
  ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
  ##
  ##========================================================
  
  distancePointLine <- function(x, y, slope, intercept) {
    ## x, y is the point to test.
    ## slope, intercept is the line to check distance.
    ##
    ## Returns distance from the line.
    ##
    ## Returns 9999 on 0 denominator conditions.
    x1 <- x-10
    x2 <- x+10
    y1 <- x1*slope+intercept
    y2 <- x2*slope+intercept
    distancePointSegment(x,y, x1,y1, x2,y2)
  }
  
  distancePointSegment <- function(px, py, x1, y1, x2, y2) {
    ## px,py is the point to test.
    ## x1,y1,x2,y2 is the line to check distance.
    ##
    ## Returns distance from the line, or if the intersecting point on the line nearest
    ## the point tested is outside the endpoints of the line, the distance to the
    ## nearest endpoint.
    ##
    ## Returns 9999 on 0 denominator conditions.
    lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
    ans <- NULL
    ix <- iy <- 0   # intersecting point
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if(any(lineMag < 0.00000001)) { # modified for vectorization by BAH
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)
    if(any((u < 0.00001) || (u > 1))) { # BAH added any to vectorize
      ## closest point does not fall within the line segment, take the shorter distance
      ## to an endpoint
      ix <- lineMagnitude(px, py, x1, y1)
      iy <- lineMagnitude(px, py, x2, y2)
      if(ix > iy)  ans <- iy
      else ans <- ix
    } else {
      ## Intersecting point is on the line, use the formula
      ix <- x1 + u * (x2 - x1)
      iy <- y1 + u * (y2 - y1)
      ans <- lineMagnitude(px, py, ix, iy)
    }
    ans
  }
  
  # End of helper functions by PB
  
  ### Now for the actual findElbow function!
  
  # Find the elbow using the method described in
  # stackoverflow.com/a/2022348/633251
  # but translated to R (see above).
  
  # Add an index to argument values for easy plotting
  DF <- data.frame(x = 1:length(y), y = y)
  fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
  m <- coef(fit)[2]
  b <- coef(fit)[1]
  
  # Check to make sure the data is concave as described
  # in the documentation, as arbitrary trends could give
  # misleading answers.  The following approach simply
  # checks to make sure all values are either above or
  # below the reference line.  This allows the values
  # to vary quite a bit and still return an answer.
  
  concave <- FALSE
  use <- 2:(nrow(DF)-1)
  refpts <- m*DF$x[use] + b
  if (all(refpts > DF$y[use]) | all(refpts < DF$y[use])) concave <- TRUE
  if (!concave) stop("Your curve doesn't appear to be concave")
  
  # Calculate the orthogonal distances
  use <- 2:(nrow(DF)-1)
  elbowd <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
  DF$dist <- c(NA, elbowd, NA) # first & last points don't have a distance
  
  if (plot) {
    edm <- which.max(DF$dist)
    plot(DF[,1:2], type = "b", xlab = "index", ylab = "y values",
         main = "Looking for the Elbow")
    segments(DF$x[1], DF$y[1],
             DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
    points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")	
    points(DF$x[edm], DF$y[edm], pch = 20)
  }
  
  if (returnIndex) return(which.max(DF$dist))
  if (!returnIndex) return(DF)
  
} # end of findElbow

elb <- PF$factors
elb$Factors <- 1:nrow(elb)
y <- elb$ExpVar
x <- elb$Factors
par(mfrow=c(1,2))
plot(x,y,type="l")
lo <- loess(y~x)
#xl <- seq(min(x),max(x), (max(x) - min(x))/19) # set the denominator to (nfactors - 1)
out = predict(lo,x)
lines(x, out, col='red', lwd=2)

findElbow(out, plot = TRUE, returnIndex = TRUE)
par(mfrow=c(1,1))


### RERUN WITHOUT COMBINING REPLICATES: 

#Importing Genus, OTU Table, Metadata, and Phylogenetic Tree: 
#Phylogeny 2: 

import_taxa <- read.table('your/directory/phyloseq-data/taxonomy.tsv', header = TRUE, sep = '\t')
head(import_taxa)

# First we have to provide names for the new columns
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxonomy <- import_taxa %>%
  mutate_at('Taxon',str_replace_all, "[a-z]__","") %>%
  separate(Taxon, sep = ';', into=ranks,remove = TRUE) %>%
  column_to_rownames(var = "Feature.ID") %>%
  as.matrix()
head(taxonomy)

TAX = tax_table(taxonomy)

#OTU Table import 2 (MAIN): 
import_table <- read.table('your/directory/phyloseq-data/feature-table.txt', header = TRUE, sep = '\t', row.names = 1, comment.char = "")
otumat <- as.matrix(import_table)
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)

#Metadata
metadata <- read.table("your/directory/DaphniaMetaData123.txt", header = T, sep = '\t')

#Read the metadata file into phyloseq:
metadata[] <- lapply(metadata, as.factor)
head(metadata)
row.names(metadata) <- paste0('X', gsub('-','.', metadata$SampleID))
tail(metadata)
META <- sample_data(metadata)

#Rooted Phylogenetic Tree:
otu_tree <- read.tree(file = "your/directory/phyloseq-data/tree.nwk")
plot(otu_tree)

#Phyloseq:
physeq <- phyloseq(OTU,TAX,META,otu_tree)
physeq <- prune_taxa(taxa_sums(physeq) > 5, physeq)
physeq


#Generate rarefaction quantiles:
qunatiles <- seq(0,1,0.05)
percentiles <- quantile(sample_sums(physeq), qunatiles)
print(percentiles)


#Rarefication Curve:
tab <- otu_table(physeq)
class(tab) <- "matrix"
tab <- t(tab)
rare <- rarecurve(tab, step = 100, ylab="ASV", label=F)



#Find what data is below a certain threshold:
sample_sums_values <- sample_sums(physeq)
threshold <- 10000
below_threshold_ids <- names(sample_sums_values[sample_sums_values < threshold])
print(below_threshold_ids)

#Rarefy:
physeq.rarefied <- rarefy_even_depth(physeq, rngseed=1, sample.size= 10000)
plot_bar(physeq.rarefied)
physeq.rarefied
sample_data(physeq.rarefied)

#Alpa Diversity Metics: 

alpha_shannon <- estimate_richness(physeq.rarefied, measures = "Shannon")

alpha_analysis <- lm(alpha_shannon$Shannon ~ Treatment, data = metaData123)
summary(alpha_analysis)

#Bray_Curtis
distance_matrix_bray <- phyloseq::distance(physeq.rarefied, method = "bray")
#treatment = shape, clone = color

#Ordination Alignment:
GP.ord <- ordinate(physeq.rarefied, "PCoA", "bray")
p1 = plot_ordination(physeq.rarefied, GP.ord, type = "taxa", color = "Genus", title = "taxa")
p1 + facet_wrap(~Phylum, 3)
print(p1)
#Metadata Alignment:
GP.ord <- ordinate(physeq.rarefied, "PCoA", "bray")
p1 = plot_ordination(physeq.rarefied, GP.ord, label = "Clone", type = "samples", color = "Treatment", title = "Treatment Exposure PCoA(Bray)")
p1 + facet_wrap(~Clone)
print(p1 + theme_bw())


#Weighted_Unifrac:

#Generate Distance Matrix; 
distance_matrix_unifrac_weighted <- phyloseq::distance(physeq.rarefied, method = "unifrac", weighted = TRUE)

#Oridnation:
unifrac.ord <- ordinate(physeq.rarefied, "PCoA", "unifrac", weighted = TRUE)
p3 = plot_ordination(physeq.rarefied, unifrac.ord, type = "taxa", color = "Genus", title = "taxa")
p3 + facet_wrap(~phylum, 3)
print(p3)
#Metadata:
unifrac.ord <- ordinate(physeq.rarefied, "PCoA", "unifrac", weighted = TRUE)
p3 = plot_ordination(physeq.rarefied, unifrac.ord, label = "Clone" , type = "samples", color = "Treatment", title = "Treatment Exposure PCoA (Weighted-Unif)")
p3 + facet_wrap(~Clone)
print(p3  + theme_bw())


#Pemranovas: #### DO ALPHA DIVERSITY
metaData123 <- sample_data(physeq.rarefied)
metaData123 <- as.matrix.data.frame(metaData123)
metaData123 <- as.data.frame(metaData123)
metaData123$NumberofGuts <- as.numeric(metaData123$NumberofGuts)
metaData123$Clone <- as.numeric(metaData123$Clone)


weighted_unifrac <- adonis2(distance_matrix_unifrac_weighted ~ Treatment + (1|Clone), data = metaData123, permutations = 1000)
print(weighted_unifrac)

bray_pr <- adonis2(distance_matrix_bray  ~ Treatment + (1|Clone), data = metaData123, permutations = 1000)
print(bray_pr)





