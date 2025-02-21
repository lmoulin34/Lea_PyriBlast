
#                              ANALYSE MICROBOME DADA2 & PHYLOSEQ (base don Callahan et al. 2016,  https://doi.org/10.1038/nmeth.3869.)

# First install a recent version of Rstudio: https://www.rstudio.com/products/rstudio/download/#download

# Install last version of R : https://cran.r-project.org/bin/windows/base/
#Select it in Tools Global Options Rversion in Rstudio.

#Installation of packages for Microbiome analyses usig DADA2

#command to install using BioManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(dada2); packageVersion("dada2")

BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("Biostrings")
BiocManager::install("ggplot2")
BiocManager::install("dplyr")
BiocManager::install("dada2")
BiocManager::install("DESeq2")
BiocManager::install("msa")
BiocManager::install("phangorn")
BiocManager::install("devtools")
BiocManager::install("Tax4Fun")
BiocManager::install("Rtools")
install.packages("Rcpp")
BiocManager::install("microbial")
#Barplots & Venn diagram packages
install.packages("microbiome")


#Load libraries

library(phyloseq)
library(DECIPHER)
library(dada2)
library(vegan)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(DESeq2)
library("microbial")


#starting the dada2 pipeline

setwd("C:/Users/moulinl/Documents/Thèse_Lea_Jobert/")


path <- "C:/Microbiote_BipolarisBlast/16S"    #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)



# Forward and reverse fastq filenames have format: SAMPLENAME-16S_1.fastq and SAMPLENAME-16S_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Find and remove primers (not compulsory, go direct to next step as we can trim them automatically by size)
# Identify primers: 337F' 805R' (size = For 21 bp, Rev 20 bp) Macrogen primers for V3V4

#FWD <- "GACTCCTACGGGAGGCWGCAG"  ## CHANGE ME to the forward primer sequence
#REV <- "GACTACCAGGGTATCTAATC"  ## CHANGE ME the reverse primer sequence



#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

#Filter and trim: Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter:
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(21, 20), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
# On Windows set multithread=FALSE; attention ici on trime les primers de 20 pb en amont de chaque sequence for et Rev (trimleft)


head(out) 
#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, pool=FALSE, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
## 251 252 253 254 255 
##   1  88 196   6   2

#Remove non-target-length sequences from your sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to "cutting a band" in-silico to get amplicons of the targeted length. 
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 379:432]


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## [1]  20 232
sum(seqtab.nochim)/sum(seqtab)
## [1] 0.964263

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "Stat_assembly16S-LeaBipoPyri.csv")



#Assign taxonomy (download first silva_nr_v138_train_set.fa.gz from silva website)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_nr_v138_train_set.fa.gz", multithread=TRUE)


#To add species information (download first silva_species_assignment_v138.fa.gz from silva website)
taxa <- addSpecies(taxa2, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_species_assignment_v138.fa.gz")


#inspect the taxonomic assignments:
taxa.print <- taxa3 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa.print <- taxa3 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#load metadata file
samdf <- read.csv2("C:/Microbiote_BipolarisBlast/metadata_tout_C1.csv")

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE pour verifier que les noms concordent

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Pathogen", "Disease", "PathogenHealth") # change with your column names in metadata samdf file
samdf <- samdf[rownames(seqtab.nochim), keep.cols]

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


# check correespondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE

#save in RDS format the seqtab  (can be reuploaded easily using redaRDS(file="") format)

saveRDS(seqtab.nochim, "16S_seqtabnochim.RDS")
saveRDS(taxa, "16S_taxa.RDS")

# extract info in a single table (without metadata)
seqtab4 <-as.data.frame(t(seqtab.nochim))
seqtab4$ID <- rownames(seqtab4)
taxa2<-as.data.frame(taxa)
taxa2$ID <- rownames(taxa)
all_data_16S  <- merge(taxa2,seqtab4,by='ID')
write.csv(all_data_16S, "all_data_16S.csv")


#phyloseq object without phylogenetic tree

ps_16S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa3))

ps_16S


#Add ASV numbering instead of fasta sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps_16S))
names(dna) <- taxa_names(ps_16S)
ps_16S  <- merge_phyloseq(ps_16S, dna)
taxa_names(ps_16S) <- paste0("ASV", seq(ntaxa(ps_16S)))
ps_16S


#remove mitochondria and chloroplast 
ps_16S <- ps_16S %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 

#remove Eukaryota
ps_16S <- ps_16S %>% subset_taxa( Kingdom!= "Eukaryota" | is.na(Kingdom) & Class!="Chloroplast" | is.na(Class) ) 
ps_16S

#alpha diversity
plot_richness(ps, x="sample_id", measures=c("Shannon", "Simpson"), color="medium")
plot_richness(ps, x="condition", measures=c("Shannon", "Simpson"), color="medium")

#Ordonates samples by condition with sortby:
plot_richness(ps_16S, x="Pathogen", measures=c("PathogenHealth", "Shannon"), color="PathogenHealth", sortby = "Shannon") + geom_boxplot()


#Remove gray font
plot_richness(ps_16S, x="Pathogen", measures=c("PathogenHealth", "Shannon"), color="PathogenHealth", sortby = "Shannon") + geom_boxplot() + theme_bw()


#Export a table with richness scores
rich = estimate_richness(ps_16S)
rich

#plot bars
plot_bar(ps_16S, fill="Phylum")


#Filter on abundance >10 reads (sup to 9 indicate it keeps 10, if you put sup 10 it means it will cut ASV at 10 and keep from 11)

ps_16S_filter10 = filter_taxa(ps_16S, function(x) sum(x > 9) > (0.0005*length(x)), TRUE)
ps_16S_filter10


#Rarefy at 10000
set.seed(1234)
ps_16S_10000F <- rarefy_even_depth(ps_16S, sample.size = 10000,
                                   rngseed = F, replace = TRUE) 
ps_16S_10000F                               

#-----------------------------------------------------
#Plot Bars (library "microbial")
library("microbial")

plotbar(ps_16S_filter10, level = "Phylum", color = NULL,group = "Pathogen",
        top = 10,return = FALSE,fontsize.x = 5,  fontsize.y = 12)


#Top 30 phylum
p1 <- plotbar(ps_16S_filter10, level = "Phylum", color = NULL,group = "Pathogen",
              top = 30,return = FALSE,fontsize.x = 5,  fontsize.y = 12)
#include a desired order f samples:
p1$data$Pathogen <- factor(p1$data$Pathogen, levels = desired_order)
print(p1)

#Top25 class
p2 <- plotbar(ps_16S_filter10, level = "Class", color = NULL,group = "Pathogen",
              top = 25,return = FALSE,fontsize.x = 10,  fontsize.y = 12)
p2$data$Pathogen <- factor(p2$data$Pathogen, levels = desired_order)
print(p2)

#Top 25 Order
p3 <- plotbar(ps_16S_filter10, level = "Order", color = NULL,group = "Pathogen",
              top = 25,return = FALSE,fontsize.x = 10,  fontsize.y = 12)
p3$data$Pathogen <- factor(p3$data$Pathogen, levels = desired_order)
print(p3)


#top25 genus
p4 <- plotbar(ps_16S_filter10, level = "Genus", color = NULL,group = "Pathogen",
              top = 25,return = FALSE,fontsize.x = 10,  fontsize.y = 12)
p4$data$Pathogen <- factor(p4$data$Pathogen, levels = desired_order)
print(p4)

#figure plot compilation
cowplot::plot_grid(p1, p2, p3, p4,nrow = 2, align = 'v', axis = 'lr', rel_heights = c(0.1, 0.1))

