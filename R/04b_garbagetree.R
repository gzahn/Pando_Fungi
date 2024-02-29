library(phyloseq)
library(Biostrings)
library(ShortRead)
library(tidyverse)
source("./R/GarbageTree_functions.r")

ps <- readRDS("./Output/clean_phyloseq_object.RDS")

seqs <- taxa_names(ps) %>% DNAStringSet()
seqnames <- paste0("ASV_",1:length(seqs))
names(seqs) <- seqnames
writeFasta(seqs,"./Output/ASV_Seqs.fasta")

# read in fasta file as named list
# file MUST be in sequential (NOT INTERLEAVED) format
seqs <- read.fastx("./Output/ASV_Seqs_oneline.fasta", type="fasta")

# cluster at multiple levels, storing result as a list of clustering tables
# 'interval' is the difference between clustering levels. i.e. 0.01 means 
# that clustering changes by 1% each time (99, 98, 97, 96...)
interval_ex <- 0.01
# startlvl is the initial (highest) clustering level
startlvl_ex <- 0.99
# finallvl is the final (lowest) clustering level
finallvl_ex <- 0.45
uc_table_list <- iterative_cluster(seqs, startlvl=startlvl_ex, finallvl=finallvl_ex, interval=interval_ex)

# transform list of clustering tables into a tree
gtree <- make_garbage_tree(uc_table_list, interval=interval_ex)

# write tree to file
write.tree(gtree, file="garbagetree_fungi.nwk")
