library(tidyverse)
library(geosphere)

source("./R/functions.R")

meta <- read_csv("./Data/pando_sample_metadata_filled.csv")
border <- read_csv("./Data/2019_border_points_latlon.csv") %>%
  dplyr::select(lon,lat)
fence <- read_csv("./Data/fence_points_latlon.csv") %>% 
  dplyr::select(lon,lat)


class(border)

pando_point <- data.frame(lon = meta$lon, lat = meta$lat)



meta <- meta %>% 
  mutate(distance_from_fence = find_gps_dists(pando_point, fence))


saveRDS(meta, "./Data/clean_Metadata.RDS")