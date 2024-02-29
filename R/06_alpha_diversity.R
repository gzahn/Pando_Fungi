# Alpha diversity

# SETUP ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(broom)
options(scipen = 999)
# functions
add_alphadiv_measures <- function(physeq){
  # build species-level physeq
  ps_species <- 
    tax_glom(physeq,"Species")
  
  alpha_asv <- estimate_richness(physeq,measures = c("Observed","Shannon","Simpson"))
  alpha_spp <- estimate_richness(ps_species,measures = c("Observed","Shannon","Simpson"))
  
  # add to sample metadata
  physeq@sam_data$asv_richness <- alpha_asv$Observed
  physeq@sam_data$asv_shannon <- alpha_asv$Shannon
  physeq@sam_data$asv_simpson <- alpha_asv$Observed
  
  # build species-level physeq
  ps_species <- 
    tax_glom(physeq,"Species")
  
  alpha_asv <- estimate_richness(physeq,measures = c("Observed","Shannon","Simpson"))
  alpha_spp <- estimate_richness(ps_species,measures = c("Observed","Shannon","Simpson"))
  
  # add to sample metadata
  physeq@sam_data$asv_richness <- alpha_asv$Observed
  physeq@sam_data$asv_shannon <- alpha_asv$Shannon
  physeq@sam_data$asv_simpson <- alpha_asv$Observed
  
  physeq@sam_data$spp_richness <- alpha_spp$Observed
  physeq@sam_data$spp_shannon <- alpha_spp$Shannon
  physeq@sam_data$spp_simpson <- alpha_spp$Observed
  
  return(physeq)
}


## data ####
# asv-level physeq
ps <- readRDS("./Output/clean_phyloseq_object.RDS")
# species-level physeq
ps_species <- 
  tax_glom(ps,"Species")


# richness
# shannon div
# at different taxonomic levels
# difference between endo/epi?

# CALCULATE ALPHA-DIV ####

# add to sample metadata
ps <- add_alphadiv_measures(ps)

# extract sample data
sam <- ps@sam_data %>% as("data.frame")

# check correlation of vars
glm(data=sam,formula = distance_from_edge ~ sample_type) %>% 
  summary

m1 <- glm(data=sam,formula = asv_richness ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "asv_richness")
m2 <- glm(data=sam,formula = asv_shannon ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "asv_shannon")
m3 <- glm(data=sam,formula = asv_simpson ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "asv_simpson")
m4 <- glm(data=sam,formula = spp_richness ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "spp_richness")
m5 <- glm(data=sam,formula = spp_shannon ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "spp_shannon")
m6 <- glm(data=sam,formula = spp_simpson ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "spp_simpson")

mod_tab <- list(m1,m2,m3,m4,m5,m6) %>% 
  reduce(full_join) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = term %>% str_remove("sample_type"))

mod_tab %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic() %>% 
  kableExtra::column_spec(5,bold = if_else(mod_tab$p.value < 0.05,TRUE,FALSE))

