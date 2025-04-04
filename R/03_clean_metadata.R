# # SETUP ####
# 
# ## Packages ####
# library(tidyverse)
# library(geosphere)
# library(janitor)
# 
# ## Functions ####
# source("./R/functions.R")
# '%ni%' <- Negate('%in%')
# ## Data ####
# meta <- read_csv("./Data/Combined_Metadata_2025.csv")
# # border <- read_csv("./Data/2019_border_points_latlon.csv") %>% 
# #   dplyr::select(lon,lat)
# # fence <- read_csv("./Data/fence_points_latlon.csv") %>% 
# #   dplyr::select(lon,lat)
# 
# # ADD EDGE DISTANCES ####
# # 2 data frames, each with only 2 cols: lon, lat (in that order)
# # pando_points <- data.frame(lon=meta$lon,
# #                            lat=meta$lat)
# # meta <- 
# # meta %>% 
# #   mutate(distance_from_edge = find_gps_dists(pando_points,
# #                                              border))
# # # ADD FENCE DISTANCES ####
# # meta <- 
# #   meta %>% 
# #   mutate(distance_from_fence = find_gps_dists(pando_points,
# #                                              fence))
# # 
# # # load manifest from seqCoast
# # manifest <- read_csv("./Data/4539_SampleManifest.csv",skip = 1) %>% clean_names()
# 
# # remove accidental duplicated sample ids
# meta2 <- meta[!duplicated(meta$sample),]
# 
# # list of file names (without path)
# # change to include pattern: "_R1_" only
# seqfiles <- list.files("./Data/Clean",pattern = "fastq.gz")
# 
# # get seqcoast id numbers for files
# l <- seqfiles %>% str_split("_")
# seqcoast_id <- map_chr(l,2)
# 
# meta2 <- 
# data.frame(seq_coast_tube_id=seqcoast_id,seqfiles) %>% 
#   mutate(seq_coast_tube_id = as.numeric(seq_coast_tube_id)) %>% 
#   full_join(manifest) %>% 
#   select(seq_coast_tube_id,seqfiles,success,sample_name) %>% 
#   rename("sample" = "sample_name") %>%  
#   full_join(meta2) %>% 
#   mutate(fwd_filepath = paste0("./Data/Raw/",seqfiles),
#          rev_filepath = paste0("./Data/Raw/",seqfiles) %>% str_replace("_R1_","_R2_")) %>% 
#   select(seq_coast_tube_id,success,sample,tree,sample_type,lat,lon,distance_from_edge,distance_from_fence,fenced,fwd_filepath,rev_filepath) %>% 
#   filter(grepl(pattern="_R1_",fwd_filepath)) %>% 
#   unique.data.frame()
# 
# meta2$sra_latlon <- 
#   paste(meta2$lat,"N",meta2$lon,"W") %>% str_remove("-")
# meta2$sra_latlon[meta2$sra_latlon %>% grep(pattern="NA")] <- NA
# meta2 <- meta2 %>% 
#   filter(sample %ni% meta2$sample[duplicated(meta2$sample)])
# 
# # SAVE CLEAN METADATA ####
# saveRDS(meta2,"./Data/Clean_Metadata.RDS")
# write_csv(meta2,"./Data/Clean_Metadata.csv")
# 
