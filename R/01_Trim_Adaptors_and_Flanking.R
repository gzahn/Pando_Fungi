# SETUP ####

## Packages ####
library(tidyverse); packageVersion('tidyverse')
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(parallel); packageVersion("parallel")

source("./R/functions.R")

## Data ####
metadata <- read_csv("./Data/Combined_Metadata_2025.csv")


# if running on local machine  
metadata$fwd_filepath <- paste0("./Data/Raw/",basename(metadata$fwd_filepath))
metadata$amplicon <- "ITS2"

## primer sequences ####

# ITS1
ITS3f_mix <- "GCATCGATGAAGAACGCAGC" 
ITS4r <- "TCCTCCGCTTATTGATATGC" 


fwd_names <- metadata$sample
list.files(file.path("./Data/Raw/filtN"))

fwd_filtn_names <- file.path("./Data/Raw/filtN",paste(fwd_names,"filtN_fwd.fastq.gz",sep="_"))
rev_filtn_names <- file.path("./Data/raw/filtN",paste(fwd_names,"filtN_rev.fastq.gz",sep="_"))
file.exists(fwd_filtn_names)

# RUN CUTADAPT ####

# PREFILTER ####

# if already done, automatically skip
filtN.fwds <- metadata$fwd_filepath[!file.exists(fwd_filtn_names)]
filtN.outs <- fwd_filtn_names[!file.exists(rev_filtn_names)]

if(length(filtN.fwds) > 0){
  filterAndTrim(fwd = filtN.fwds,
                filt = fwd_filtn_names,
                maxN = 0,
                compress = TRUE,
                multithread = parallel::detectCores())
} else {
  message("All files previously had Ns removed...skipping filtN step.")
}


## On ITS Samples ####
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(ITS3f_mix)
REV.orients <- allOrients(ITS4r)
FWD.orients; REV.orients

# Discover primer matches, regardless of orientation ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtN.outs[[1]]))
system2("cutadapt", args = "--version")

path.cut <- file.path("./Data/Raw/cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(filtN.outs))

FWD.RC <- dada2:::rc(ITS3f_mix)

cutadapt_out <- file.path(path.cut,filtN.outs %>% basename %>% str_replace("_filtN_","_cutadapt_"))
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", ITS3f_mix, "-a", FWD.RC) 

for(i in seq_along(filtN.outs)) {
  system2("cutadapt", args = c(R1.flags, "-n", 2, "--minimum-length 100", # -n 2 required to remove FWD and REV from reads
                               "-o", cutadapt_out[i], 
                               filtN.outs[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = cutadapt_out[[1]]))


# move cutadapt files to "Clean" directory
file.rename(cutadapt_out,file.path("./Data/Clean",basename(cutadapt_out)))



# Run itsxpress
  #on just ITS samples (and just FWD ITS1 samples, actually)

# run_itsxpress(directory="./data/raw/cutadapt", # where cutadapted reads live
#               itsregion="ITS1", # must be "ITS1" or "ITS2"
#               taxa_group="All",
#               nthreads=(parallel::detectCores()-1),
#               fwd_pattern="ITS_cutadapt_fwd.fastq.gz",
#               rev_pattern="ITS_cutadapt_rev.fastq.gz",
#               itsxpress.path="/home/gzahn/.local/bin/itsxpress", #path to executable
#               fwd.only=TRUE)

# 
# metadata$local_fwd_cutadapt_paths <- cutadapt_ftest_fwd
# metadata$local_rev_cutadapt_paths <- cutadapt_ftest_rev
# 
# saveRDS(metadata,"./data/cutadapt_metadata.RDS")
