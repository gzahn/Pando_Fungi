# -----------------------------------------------------------------------------#
# Pando foliar fungi ITS DADA2 pipeline
# Processing raw amplicon reads (these have already been passed through itsxpress)
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     dada2 v 1.24.0
#                     decontam v 1.16.0
#                     phyloseq v 1.40.0
#                     Biostrings v 2.64.0
#                     patchwork v 1.1.1
#                     readxl v 1.4.1
#                     janitor::clean_names() v 2.1.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal  # 
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
source("./R/functions.R")

set.seed(666) # "random" seed for reproducibility

'%ni%' <- Negate('%in%')

# PARSE FILE PATHS ####

# File parsing - 

path <- "./Data/Clean" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# change filenames since itsxpress output "ITS1"
new.filenames <- list.files(path,full.names = TRUE,pattern=".fastq.gz")
 
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "cutadapt_fwd.fastq.gz")) # make pattern match your FWD reads

sample.ids <- basename(fns) %>%  str_remove("_cutadapt_fwd.fastq.gz")
meta <- read_csv("./Data/Combined_Metadata_2025.csv")
row.names(meta) <- as.character(meta$sample)
sample.names <- meta[sample.ids,"sample"]

# some files may not have corresponding metadata for some reason!?
goodsamples <- which(!is.na(sample.names))
meta <- meta[sample.ids,][goodsamples,]
sample.ids <- sample.ids[goodsamples]
fns <- fns[goodsamples]



# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[11]) + ggtitle("Example forward reads")

# display and save the plots
p1
ggsave("./Output/Figs/unfiltered_quality_plots.png",dpi=500,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered", paste0(sample.ids, "_FWD_filt.fastq.gz"))

# this is the actual qualit control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=3, # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     compress=TRUE, # compress output files with gzip
                     multithread=parallel::detectCores()-1) # On Windows set multithread=FALSE

# if Some input samples had no reads pass the filter....
# ... clean up metadata to remaining samples.

remaining_samples <- which(out[,2] >= 1000)
meta <- meta[remaining_samples,]
write_csv(meta,"./Data/Combined_Metadata_2025_remaining_samples.csv")

# save the filtration info in case we need it later
saveRDS(out, "./Output/trackreads.RDS")

list.files(filtpath, full.names = TRUE)
# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")
sample.names <- sample.names[sample.names %in% meta$sample]

# remove cutadapt files that didn't pass threshold (>=1000 reads)
missed_cutoff <- filts_f %>% str_remove("./Data/Clean/filtered/") %>% str_remove("_FWD_filt.fastq.gz") %ni% sample.names
file.remove(filts_f[missed_cutoff])
filts_f <- filts_f[!missed_cutoff]


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[1]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[1])+ ggtitle("Filtered")
p3 / p4
ggsave("./Output/Figs/filtered_quality_comparison.png",dpi=300,height = 6,width = 6)




# LEARN ERROR RATES ####

# learn errors
errF <- learnErrors(filts_f, multithread=parallel::detectCores()-1, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errF,"./Output/errF.RDS")
errF <- readRDS("./Output/errF.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./Output/figs/error_model.png",dpi=200,height = 6,width = 6)

# dereplication
derepF <- derepFastq(filts_f, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepF) <- 
  names(derepF) %>% 
  str_remove("_FWD_filt.fastq.gz")
saveRDS(derepF,"./Output/derepF.RDS")
# derepF <- readRDS("./Output/derepF.RDS")

# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = FALSE, verbose=TRUE)
saveRDS(dadaFs,"Output/dadaFs.RDS") 

# if dada worked and was saved properly, remove giant derep
if(file.exists("Output/dadaFs.RDS") & (length(derepF) > 1)){
  rm(derepF)
  file.remove("./Output/derepF.RDS")
}
  
# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(dadaFs)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./Output/seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF","nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./Output/ITS_read_counts_at_each_step.csv", row.names = TRUE)


# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./Output/seqtab.nochim.clean.RDS")
dim(seqtab.nochim)

# filter metadata to match remaining samples
meta <- meta[meta$sample %in% row.names(seqtab.nochim),]

# Find and remove contaminants ####
contams = isContaminant(seqtab.nochim, neg = grepl("blank",meta$sample), normalize = TRUE)
table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./Output/likely_contaminants.csv", row.names = TRUE)

# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[!grepl("blank",row.names(seqtab.nochim)),]
meta = meta[!grepl("blank",meta$sample),]
seqtab.nochim <- seqtab.nochim[row.names(seqtab.nochim) %in% meta$sample,]
meta <- meta[meta$sample %in% row.names(seqtab.nochim),]
dim(seqtab.nochim)
dim(meta)

# reload point
saveRDS(seqtab.nochim,"./Output/seqtab_for_taxonomy.RDS")
seqtab.nochim <- readRDS("./Output/seqtab_for_taxonomy.RDS")

# ASSIGN TAXONOMY ####
tax <- assignTaxonomy(seqtab.nochim,
                      refFasta = "./Taxonomy/Eukaryome_General_ITS_v1.8_reformatted.fasta.gz",
                      multithread = TRUE,
                      tryRC = FALSE,
                      verbose = TRUE)
tax
beepr::beep(sound=2)
saveRDS(tax,"./Output/sequence_taxonomy.RDS")
tax <- readRDS("./Output/sequence_taxonomy.RDS")
dim(seqtab.nochim)
dim(meta)
dim(tax)

# BUILD PHYLOSEQ ####

met <- sample_data(meta)
sample_names(met) <- meta$sample
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
sample_names(otu)
taxa <- tax_table(tax)
taxa_names(taxa)
sample_names(taxa)

ps <- phyloseq(met,otu,taxa)

saveRDS(ps,"./Output/phyloseq_object_not-cleaned.RDS")

# clean up physeq
ps <- 
  pd %>% 
  clean_ps_taxonomy() %>% 
  subset_taxa(Kingdom == "Fungi")
ps <- 
  ps %>% 
  subset_taxa(taxa_sums(ps) > 1)
ps <- 
  ps %>% 
  subset_samples(sample_sums(ps) >= 100)
saveRDS(ps,"./Output/clean_phyloseq_object.RDS")
