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


# SAVE CLEAN METADATA ####
saveRDS(meta,"./Data/Clean_Metadata.RDS")

