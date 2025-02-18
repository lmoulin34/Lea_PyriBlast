# Script rarefaction curve and alpha diversity
# Project: Pathogen Impact on Root Microbiome

# Install required packages
library(dada2)
library(phyloseq)
library(dplyr)
library(data.table)
library(ggplot2)
library(vegan)
library(phangorn)
library(DECIPHER)
library(RColorBrewer)
library(decontam)
library(scales)

# Set working directory and define path
setwd("path/to/directory")
path <- "path/to/directory"
list.files(path)

# Import necessary data files
taxa <- read.csv(file = "path/to/directory/taxa.csv", header = TRUE, row.names = 1, sep = ";")
metadata <- read.csv2("path/to/directory/metadata.csv")
asvITS <- readRDS(file = "path/to/directory/seqtab.nochim.RDS")

# Ensure sample names match
all(rownames(asvITS) %in% metadata$SampleID)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[rownames(asvITS), ]

# Create phyloseq object
otu_tab <- otu_table(as.matrix(asvITS), taxa_are_rows = FALSE)
taxa_tab <- tax_table(as.matrix(taxa))
psITS <- phyloseq(otu_tab, sample_data(metadata), taxa_tab)


# Rename ASVs
dna.ITS <- Biostrings::DNAStringSet(taxa_names(psITS))
names(dna.ITS) <- taxa_names(psITS)
psITS <- merge_phyloseq(psITS, dna.ITS)
taxa_names(psITS) <- paste0("ASV", seq(ntaxa(psITS)))

# Rarefaction curve
RC <- ggrare(psITS, step = 100, se = FALSE, color = "Condition") + 
  labs(title = "Rarefaction Curve", x = "Number of ASVs", y = "Depth") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
print(RC)

# Taxonomic filtering
psITS <- subset_taxa(psITS, !is.na(Kingdom) & !is.na(Phylum))
psITS <- subset_taxa(psITS, !Kingdom %in% c("k__Viridiplantae"))

# Abundance filtering
psITS <- filter_taxa(psITS, function(x) sum(x >= 9) > 0, TRUE)
psITS <- prune_samples(sample_sums(psITS) >= 15000, psITS)


# Phylogenetic tree
seqs.ITS <- refseq(psITS)
alignment.ITS <- AlignTranslation(seqs.ITS, sense = "+", readingFrame = 2, type = "DNAStringSet")
phang.align.ITS <- phyDat(as(alignment.ITS, "matrix"), type = "DNA")
dm.ITS <- dist.ml(phang.align.ITS)
treeNJ.ITS <- NJ(dm.ITS)
treeITS <- phyloseq(tax_table(psITS), sample_data(psITS), otu_table(psITS, taxa_are_rows = FALSE), phy_tree(treeNJ.ITS))
saveRDS(treeITS, "path/to/directory/treeITS.rds")




#Table of alpha-diversity estimators
table_ITS <- estimate_richness(psITS, split = TRUE, measures=c("Observed", "Shannon"))
#Add evenness
H.ITS <- table_ITS$Shannon
S1.ITS <- table_ITS$Observed
S.ITS <- log(S1.ITS)
eve_ITS <- H.ITS/S.ITS
table_ITS$Evenness <- eve_ITS
#Bind sampledata + table of alpha_div (here R1 only)
datar.ITS <- cbind(sample_data(psITS), table_ITS)
#Add phylogenetic diversity (Faith's index)
PD.ITS <- as.data.frame.matrix(otu_table(psITS))
Faith.ITS <- pd(samp = PD.ITS, tree = phy_tree(treeITS), include.root = F) 
d.faith.ITS <- merge(sample_data(treeITS), Faith.ITS, by="row.names", all=TRUE) 



# Define your desired order
desired_order <- c("BipolarisHealthy", "BipolarisDisease", "PyriculariaHealthy", "PyriculariaDisease")
#Richness
p.ric.ITS <- ggplot(data = datar.ITS, aes(x = factor(PathogenHealth, level = desired_order), y = Observed)) + 
  geom_boxplot(outlier.shape = NA, fill = c("#0E9594", "#FC6471", "#9BC59D", "#DE9E36"), color = "black", alpha = 0.8) +
  geom_jitter(width = 0.15) +
  xlab('Condition') +
  ylab('Observed Richness') +
  ggtitle("Observed ASV Richness in Healthy or Diseased Root for Both Pathosystems (ITS)") +
  theme_minimal(base_size = 15) +  
  theme(
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12, color = "black"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",  # No legend needed
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  scale_x_discrete(labels = c("BipolarisHealthy" = "Bipolaris\nHealthy",
                              "BipolarisDisease" = "Bipolaris\nDisease",
                              "PyriculariaHealthy" = "Pyricularia\nHealthy",
                              "PyriculariaDisease" = "Pyricularia\nDisease")) +
  stat_compare_means(comparisons = list(c("BipolarisHealthy", "BipolarisDisease"),
                                        c("PyriculariaHealthy", "PyriculariaDisease")),
                     method = "wilcox.test", label = "p.signif")

# Display and save the plot
print(p.ric.ITS)
ggsave("p.ric.ITS.svg", plot = p.ric.ITS, width = 12, height = 6, device = "svg")


#Evenness
p.eve.ITS <- ggplot(data = datar.ITS, aes(x = factor(PathogenHealth, level = desired_order), y = Evenness)) + 
  geom_boxplot(outlier.shape = NA, fill = c("#0E9594", "#FC6471", "#9BC59D", "#DE9E36"), color = "black", alpha = 0.8) +
  geom_jitter(width = 0.15) +
  xlab('Condition') +
  ylab('Pielou Evenness') +
  ggtitle("Pielou Evenness in Healthy or Diseased Root for Both Pathosystems (ITS)") +
  theme_minimal(base_size = 15) +  
  theme(
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12, color = "black"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",  # No legend needed
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  scale_x_discrete(labels = c("BipolarisHealthy" = "Bipolaris\nHealthy",
                              "BipolarisDisease" = "Bipolaris\nDisease",
                              "PyriculariaHealthy" = "Pyricularia\nHealthy",
                              "PyriculariaDisease" = "Pyricularia\nDisease")) +
  stat_compare_means(comparisons = list(c("BipolarisHealthy", "BipolarisDisease"),
                                        c("PyriculariaHealthy", "PyriculariaDisease")),
                     method = "wilcox.test", label = "p.signif")

# Display and save the plot
print(p.eve.ITS)
ggsave("p.eve.ITS.svg", plot = p.eve.ITS, width = 12, height = 6, device = "svg")


#Phylogenetic Diversity

# Desired order of factor levels
desired_order <- c("BipolarisHealthy", "BipolarisDisease", "PyriculariaHealthy", "PyriculariaDisease")

# Create the enhanced plot
p.pd.ITS <- ggplot(data = d.faith.ITS, aes(x = factor(PathogenHealth, level = desired_order), y = PD)) + 
  geom_boxplot(outlier.shape = NA, fill = c("#0E9594", "#FC6471", "#9BC59D", "#DE9E36"), color = "black", alpha = 0.8) +
  geom_jitter(width = 0.15) +
  xlab('Condition') +
  ylab('Phylogenetic Diversity') +
  ggtitle("Phylogenetic Diversity in Healthy or Diseased Root for Both Pathosystems (ITS)") +
  theme_minimal(base_size = 15) +  
  theme(
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12, color = "black"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",  # No legend needed
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  scale_x_discrete(labels = c("BipolarisHealthy" = "Bipolaris\nHealthy",
                              "BipolarisDisease" = "Bipolaris\nDisease",
                              "PyriculariaHealthy" = "Pyricularia\nHealthy",
                              "PyriculariaDisease" = "Pyricularia\nDisease")) +
  stat_compare_means(comparisons = list(c("BipolarisHealthy", "BipolarisDisease"),
                                        c("PyriculariaHealthy", "PyriculariaDisease")),
                     method = "wilcox.test", label = "p.signif")

# Display and save the plot
print(p.pd.ITS)
ggsave("p.pd.ITS.svg", plot = p.pd.ITS, width = 12, height = 6, device = "svg")



