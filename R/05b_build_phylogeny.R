# -----------------------------------------------------------------------------#
# modified from Meta-amplicon analysis recipe
# Building and adding a phylogeny to the cleaned phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     vegan v 2.5.6
#                     phangorn v 2.2.5
#                     phyloseq v 1.30.0
#                     msa v 1.18.0
#                     ape v 5.4
#                     seqinr v 3.6.1
<<<<<<< HEAD
#                     phylogram v 2.1.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Perform multiple sequence alignment of all ASVs, build distance matrix,       # 
# construct and refine a phylogenetic tree, add the tree to the phyloseq object #
#           With larger data sets, this can be a long process...                #
# Further, proper phylogenetics is beyond the scope of this tutorial.           #
#################################################################################

# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
<<<<<<< HEAD
library(DECIPHER); packageVersion("DECIPHER")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(phylogram); packageVersion('phylogram')
library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")

# Read in phyloseq object from first script output ####
ps <- readRDS("./Output/clean_phyloseq_object.RDS")
sample_data(ps)

# simplify ASV names
seqs <- rownames(tax_table(ps))
<<<<<<< HEAD
# names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree
names(seqs) <- seqs
# Multiple sequence alignment  ####
# alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)
alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs))
# save progress 
saveRDS(alignment,"./Output/ITS_dna_alignment_muscle.RDS")
alignment <- readRDS("./Output/ITS_dna_alignment_muscle.RDS")
# # Convert to phangorn format
# phang.align = as.phyDat(as.character(alignment), type = "DNA")

# distance - maximum likelihood ####
# dm <- dist.ml(phang.align)
dist <- DistanceMatrix(alignment,type = "dist",correction="Jukes-Cantor", verbose=FALSE,processors = NULL)
# distm <- DistanceMatrix(alignment,type = "matrix",correction="Jukes-Cantor", verbose=FALSE,processors = NULL)
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# Multiple sequence alignment  ####
alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)

# save progress 
saveRDS(alignment,"./Output/ITS_dna_alignment_muscle.RDS")

# Convert to phangorn format
phang.align = as.phyDat(alignment, type = "DNA")

# distance - maximum likelihood ####
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"./Output/ITS_ML_Distance.RDS")

# Initial neighbor-joining tree ####
<<<<<<< HEAD
# treeNJ_phangorn <- njs(dm) # Note, tip order != sequence order
treeNJ <- TreeLine(myDistMatrix=dist, method="NJ", cutoff=0.05, showPlot=FALSE, verbose=FALSE)
# treeNJm <- TreeLine(myDistMatrix=distm, method="NJ", cutoff=0.05, showPlot=FALSE, verbose=FALSE)

# save progress
saveRDS(treeNJ, "./Output/ITS_treeNJ.RDS")
treeNJ <- readRDS("./Output/ITS_treeNJ.RDS")

phylogram::as.phylo(treeNJ)
tree <- phylogram::as.phylo(treeNJ)
tree$tip.label

taxa_names(ps)
treeNJ <- NJ(dm) # Note, tip order != sequence order

# save progress
saveRDS(treeNJ, "./Output/ITS_treeNJ.RDS")

# Estimate model parameters ####
fit = pml(treeNJ, data=phang.align)

#save
saveRDS(fit,"./Output/ITS_fit_treeNJ.RDS")


fit$tree$tip.label <- seqs

# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
<<<<<<< HEAD
                phy_tree(tree))


                phy_tree(fit$tree))


# Save updated phyloseq object with tree
saveRDS(ps2, "./Output/ps_cleaned_w_tree.RDS")

