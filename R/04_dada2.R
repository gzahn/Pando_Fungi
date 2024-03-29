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

set.seed(666) # "random" seed for reproducibility


# PARSE FILE PATHS ####

# File parsing - 

<<<<<<< HEAD
path <- "./Data/Raw" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
path <- "./Data/Clean" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# change filenames since itsxpress output "ITS1"
old.filenames <- list.files(path,full.names = TRUE,pattern=".fastq.gz")
new.filenames <- old.filenames %>% str_replace("ITS1","ITS2")
file.rename(old.filenames,new.filenames)

<<<<<<< HEAD
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_001.fastq.gz")) # make pattern match your FWD reads
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "ITS2.fastq.gz")) # make pattern match your FWD reads
sample.ids <- basename(fns) %>% str_split("_") %>% map_chr(2)
meta <- readRDS("./Data/Clean_Metadata.RDS")
row.names(meta) <- as.character(meta$seq_coast_tube_id)
sample.names <- meta[sample.ids,"sample"]

# some files may not have corresponding metadata for some reason!?
goodsamples <- which(!is.na(sample.names))
meta <- meta[sample.ids,][goodsamples,]
sample.ids <- sample.ids[goodsamples]
sample.names <- sample.names[goodsamples]
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
<<<<<<< HEAD
                     maxEE=3, # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     compress=TRUE, # compress output files with gzip
                     multithread=20) # On Windows set multithread=FALSE
                     maxEE=2, # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     compress=TRUE, # compress output files with gzip
                     multithread=16) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./Output/trackreads.RDS")


# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[1]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[1])+ ggtitle("Filtered")
p3 / p4
ggsave("./Output/Figs/filtered_quality_comparison.png",dpi=300,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
<<<<<<< HEAD
errF <- learnErrors(filts_f, multithread=16, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errF,"./Output/errF.RDS")
errF <- readRDS("./Output/errF.RDS")
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2) # set multithread = FALSE on Windows
saveRDS(errF,"./Output/errF.RDS")

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
<<<<<<< HEAD
saveRDS(derepF,"./Output/derepF.RDS")
derepF <- readRDS("./Output/derepF.RDS")

# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = FALSE, verbose=TRUE) # set multithread = FALSE on Windows
saveRDS(dadaFs,"Output/dadaFs.RDS")
rm(derepF)
# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(dadaFs)
rm(dadaFs)

# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
saveRDS(dadaFs,"Output/dadaFs.RDS")


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
meta <- meta[row.names(meta) %in% row.names(seqtab.nochim),]

# Find and remove contaminants ####
contams = isContaminant(seqtab.nochim, neg = grepl("blank",meta$sample), normalize = TRUE)
table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./Output/likely_contaminants.csv", row.names = TRUE)

# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[!grepl("blank",meta$sample),]
meta = meta[!grepl("blank",meta$sample),]
dim(seqtab.nochim)


# reload point
saveRDS(seqtab.nochim,"./Output/seqtab_for_taxonomy.RDS")

# ASSIGN TAXONOMY ####
tax <- assignTaxonomy(seqtab.nochim,
                      refFasta = "./Taxonomy/sh_general_release_dynamic_all_18.07.2023_dev.fasta.gz",
                      multithread = TRUE,
                      tryRC = TRUE,
                      verbose = TRUE)
tax
beepr::beep(sound=2)
<<<<<<< HEAD
saveRDS(tax,"./Output/sequence_taxonomy.RDS")

dim(seqtab.nochim)
dim(meta)
dim(tax)

# BUILD PHYLOSEQ ####

met <- sample_data(meta)
sample_names(met)
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
sample_names(otu)
taxa <- tax_table(tax)
<<<<<<< HEAD
taxa_names(taxa)
sample_names(taxa)

ps <- phyloseq(met,otu,taxa)

saveRDS(ps,"./Output/phyloseq_object_not-cleaned.RDS")
<<<<<<< HEAD
ps

