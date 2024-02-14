# SETUP ####

## Packages ####
library(tidyverse)
library(geosphere)

## Functions ####
source("./R/functions.R")

## Data ####
meta <- read_csv("./Data/pando_sample_metadata_filled.csv")
border <- read_csv("./Data/2019_border_points_latlon.csv") %>% 
  dplyr::select(lon,lat)
fence <- read_csv("./Data/fence_points_latlon.csv") %>% 
  dplyr::select(lon,lat)

# ADD EDGE DISTANCES ####
# 2 data frames, each with only 2 cols: lon, lat (in that order)
pando_points <- data.frame(lon=meta$lon,
                           lat=meta$lat)
meta <- 
meta %>% 
  mutate(distance_from_edge = find_gps_dists(pando_points,
                                             border))
# ADD FENCE DISTANCES ####
meta <- 
  meta %>% 
  mutate(distance_from_fence = find_gps_dists(pando_points,
                                             fence))

# load manifest from seqCoast
manifest <- read_csv("./Data/4539_SampleManifest.csv",skip = 1) %>% clean_names()

# remove accidental duplicated sample ids
meta <- meta[!duplicated(meta$sample),]

# list of file names (without path)
# change to include pattern: "_R1_" only
seqfiles <- list.files("./Data/Seqs/")

# get seqcoast id numbers for files
l <- seqfiles %>% str_split("_")
seqcoast_id <- map_chr(l,2)

meta <- 
data.frame(seq_coast_tube_id=seqcoast_id,seqfiles) %>% 
  mutate(seq_coast_tube_id = as.numeric(seq_coast_tube_id)) %>% 
  full_join(manifest) %>% 
  select(seq_coast_tube_id,seqfiles,success,sample_name) %>% 
  rename("sample" = "sample_name") %>%  
  full_join(meta)

meta <- meta[!duplicated(meta$sample),]


# SAVE CLEAN METADATA ####
saveRDS(meta,"./Data/Clean_Metadata.RDS")

