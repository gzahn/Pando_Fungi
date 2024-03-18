# SETUP ####
library(tidyverse)
library(phyloseq)
library(patchwork)
library(vegan)
library(Biostrings)

set.seed(666)
'%ni%' <- Negate('%in%')
ra <- function(x){x/sum(x)}

# LOAD DATA ####
ps <- readRDS("./Output/phyloseq_object_not-cleaned.RDS")
ps
# CLEAN PS OBJECT ####


# Subset to fungi only
fung <- subset_taxa(ps,Kingdom == "k__Fungi")




# remove low-count samples
fung <- subset_samples(fung, sample_sums(fung) >= 10)

# remove low-count ASVs
taxa_sums(fung) %>% plot
fung <- subset_taxa(fung, taxa_sums(fung) > 1)

# not sure why these don't work
# add refseq slot
# ps@refseq <- taxa_names(ps) %>% DNAStringSet()

# clean up taxonomy assignments
# for(i in 1:length(rank_names(ps))){
#   ps@tax_table[,i] <- ps@tax_table[,i] %>% str_remove(".__")
##   names(ps@tax_table[,i]) <- ps@refseq
# }


# remove low-count samples
fung <- subset_samples(fung, sample_sums(fung) >= 100)

# remove low-count ASVs
taxa_sums(fung) %>% plot
fung <- subset_taxa(fung, taxa_sums(fung) > 0)


# remove taxa with no phylum assignment
fung <- subset_taxa(fung, !is.na(Phylum))

# EXPORT ####
saveRDS(fung,"./Output/clean_phyloseq_object.RDS")


##################################################
fung %>% 
  merge_samples("sample_type") %>% 
  transform_sample_counts(ra) %>% 
  plot_bar(fill="Class")
ps@sam_data %>% head



fung %>% 
  transform_sample_counts(ra) %>% 
  plot_bar(fill="Phylum")

# questions...
# is endophyte community subset of epiphytes? (nestedness tests)
# is there a spatial component to this pattern?
# different inside/outside of fence?
# quantifying ordination / beta-diversity

ord <- 
fung %>% 
  transform_sample_counts(ra) %>% 
  ordinate(method = "NMDS")

plot_ordination(fung, ord, color="distance_from_edge") +
  facet_wrap(~sample_type,scales = 'free')

