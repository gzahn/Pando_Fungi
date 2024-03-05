# setup ####
library(phyloseq)
library(tidyverse)

'%ni%' <- Negate('%in%')

# load data ####
ps <- readRDS("./Output/ps_cleaned_w_tree.RDS")

ps@sam_data$sample_type %>% table
# relax parameters in cleaning step to avoid losing so many samples???

# for each tree where there are both epi and endo...
# calculate proportion of endo taxa also found in epi sample

prop_nested <- c()
for(i in ps@sam_data$tree){
  x <- ps %>% subset_samples(tree == i)
  x <- x %>% subset_taxa(taxa_sums(x) > 0)
  if(nsamples(x) < 2){prop_nested[i] <- NA;next}
  endo <- x %>% subset_samples(sample_type == "endophyte")
  endo <- endo %>% subset_taxa(taxa_sums(endo) > 0)
  epi <- x %>% subset_samples(sample_type == "epiphyte")
  epi <- epi %>% subset_taxa(taxa_sums(epi) > 0)
  prop_nested[i] <- sum(taxa_names(endo) %in% taxa_names(epi)) / ntaxa(endo)
}
prop_nested %>% summary
