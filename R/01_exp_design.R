library(tidyverse)
library(edibble)
library(readxl)
# Build experimental design

samples <- read_xlsx("./Data/pando_pcr.xlsx",sheet = 2)

n.samples <- samples[complete.cases(samples),] %>% nrow()

exp <- 
  design("Pando endophytes and epiphytes")
exp <- 
exp %>% 
  set_units(tree=96,
            sample = nested_in(tree,2)
            ) %>% 
  set_rcrds_of(tree = c("lat","lon","fenced","distance_from_edge"),
               sample = c('sample_type','plate','row','col','fwd_filepath','rev_filepath')
               ) %>% 
  expect_rcrds(factor(sample_type, levels = c("Endophyte","Epiphyte")),
               factor(fenced, levels = c("TRUE","FALSE")),
               distance_from_edge >= 0,
               )
tree_id <- 
samples$sample %>% 
  str_to_lower() %>% 
  str_remove("end_") %>% str_remove("epi_") %>% 
  str_replace("blank_1","NA") %>% str_replace("blank_2","NA") %>% str_replace("blank1","NA")

sample.type <- 
  samples$sample %>% str_split("_") %>% map_chr(1)

tab <- 
exp %>% 
  edibble::serve_table() %>% 
  mutate(sample = samples$sample %>% str_to_lower(),
         plate = samples$plate,
         row = samples$row,
         col = samples$col,
         tree = tree_id,
         sample_type = sample.type,
         sample_type = case_when(sample_type == "epi" ~ "epiphyte",
                                 sample_type == "end" ~ "endophyte"))
  
tab
plot(exp)
tab %>% 
  write_csv("./Data/pando_sample_metadata.csv")

