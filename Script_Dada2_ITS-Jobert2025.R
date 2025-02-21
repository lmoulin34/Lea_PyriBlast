
#                              ANALYSE MICROBOME DADA2 & PHYLOSEQ with ITS



#Load libraries (install packages first if needed)

library(phyloseq)
library(ShortRead)
library(DECIPHER)
library(dada2)
library(vegan)
library(Biostrings)
library(ggplot2)
library(dplyr)



#starting the dada2 pipeline
setwd("C:/Microbiote_BipolarisBlast/ITS")

path <- "C:/Microbiote_BipolarisBlast/ITS"    #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


#DETECT AND REMOVE PRIMERS
# Forward and reverse fastq filenames have format: SAMPLENAME-ITS_1.fastq and SAMPLENAME-ITS_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names


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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(20, 20), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
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
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 260:433]


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
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

write.csv(track, "Stat_assemblyDADA2ITS.csv")


#Assign taxonomy (download first silva_nr_v132_train_set.fa.gz from silva website)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Scripts R/Taxonomic_databases/sh_general_release_dynamic_s_all_25.07.2023_dev.fasta", multithread=TRUE)


#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

View(sample.names)

#load metadata file
samdf <- read.csv2("C:/Microbiote_BipolarisBlast/metadata.csv")

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE to verify names are the sames in samdf and seqtab file

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Inoc") # give names of metadata columns in samdf you want to keep 
samdf <- samdf[rownames(seqtab.nochim), keep.cols]

#samples.out <- rownames(seqtab.nochim)
#rownames(samdf) <- samples.out


# check correespondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE

#save in RDS format the seqtab - format you can reupload later with readRDS(file="")  function

saveRDS(seqtab.nochim, "seqtab.nochim.RDS")
saveRDS(taxa, "taxa_Lea_ITSBiPy.RDS")


# extract info in a single table (without metadata)
seqtab3 <-as.data.frame(t(seqtab.nochim))
seqtab3$ID <- rownames(seqtab3)
taxa2<-as.data.frame(taxa)
taxa2$ID <- rownames(taxa)
all_data  <- merge(taxa2,seqtab3,by='ID')
write.csv(all_data, "all_data_ITS.csv")  #Export table in csv format, this table can be open in excel to work 

#phyloseq object without phylogenetic tree

ps_its <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_its


