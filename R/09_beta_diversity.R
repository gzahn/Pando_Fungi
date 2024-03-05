# SETUP ####
library(phyloseq)
library(vegan)
library(tidyverse)

# data ####
ps <- readRDS("./Output/ps_cleaned_w_tree.RDS")
ps_ra <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)})
ord <- ordinate(ps_ra,method = "NMDS",distance = 'unifrac')
plot_ordination(ps,ord,color="sample_type") + stat_ellipse()
plot_ordination(ps,ord,color="distance_from_edge")

# Permutational ANOVA
adonis2(ps_ra@otu_table ~ ps_ra@sam_data$sample_type + ps_ra@sam_data$distance_from_edge)

# OUTCOME ~ PREDICTORS