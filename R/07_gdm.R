library(phyloseq)
library(tidyverse)
library(gdm)


fung <- readRDS("./Output/clean_phyloseq_object.RDS")
# extract species by site info
epi <- 
  fung %>% 
  subset_samples(sample_type == "epiphyte")
epi <- 
  epi %>% 
  subset_taxa(taxa_sums(epi) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)})
epi <- psmelt(epi)
endo <- 
  fung %>% 
  subset_samples(sample_type == "endophyte")
endo <- 
  endo %>% 
  subset_taxa(taxa_sums(endo) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)})
endo <- psmelt(endo)
# biological data
# get columns with xy, site ID, and species data
sppTab_epi <- epi %>% select(OTU,tree,lon,lat,Abundance)
sppTab_endo <- endo %>% select(OTU,tree,lon,lat,Abundance)
# get columns with site ID, env. data, and xy-coordinates
envTab_epi <- epi %>% select(tree, lon, lat, distance_from_edge, distance_from_fence)
envTab_endo <- endo %>% select(tree, lon, lat, distance_from_edge, distance_from_fence)
# format for gdm
gdmTab_epi <- formatsitepair(bioData=sppTab_epi, 
                             bioFormat=2, #x-y spp list
                             XColumn="lon", 
                             YColumn="lat",
                             sppColumn="OTU", 
                             siteColumn="tree", 
                             predData=envTab_epi,
                             abundance = TRUE,
                             abundColumn = "Abundance")
gdmTab_endo <- formatsitepair(bioData=sppTab_endo, 
                              bioFormat=2, #x-y spp list
                              XColumn="lon", 
                              YColumn="lat",
                              sppColumn="OTU", 
                              siteColumn="tree", 
                              predData=envTab_endo,
                              abundance = TRUE,
                              abundColumn = "Abundance")

# fit GDM models
gdm_epi <- gdm(data = gdmTab_epi,geo = TRUE)
gdm_endo <- gdm(data = gdmTab_endo,geo = TRUE)
# quick look at model fits
summary(gdm_epi)
summary(gdm_endo)

plot(gdm_epi)
plot(gdm_endo)

# predictions from model (using same distances)
gdm_epi_pred <- predict(object=gdm_epi, data=gdmTab_epi)
gdm_endo_pred <- predict(object=gdm_endo, data=gdmTab_endo)

epi_preds <- data.frame(observed = gdmTab_epi$distance,
                        predicted = gdm_epi_pred,
                        sample_type = "epiphyte")
endo_preds <- data.frame(observed = gdmTab_endo$distance,
                         predicted = gdm_endo_pred,
                         sample_type = "endophyte")
full_join(epi_preds,endo_preds) %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~sample_type,scales = 'free')
plot(gdm_endo)
plot(gdm_epi)
