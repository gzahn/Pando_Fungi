library(phyloseq)
library(tidyverse)
library(gdm)
library(patchwork)
set.seed(666)
sampletypecolors <- c("darkblue","#fca311")

fung <- readRDS("./Output/ps_cleaned_w_tree.RDS")
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
envTab_epi <- epi %>% select(tree, lon, lat, distance_from_edge)
envTab_endo <- endo %>% select(tree, lon, lat, distance_from_edge)
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
                        dist_from_edge = gdm_epi$ecological,
                        sample_type = "epiphyte")
endo_preds <- data.frame(observed = gdmTab_endo$distance,
                         predicted = gdm_endo_pred,
                         dist_from_edge = gdm_endo$ecological,
                         sample_type = "endophyte")
full_join(epi_preds,endo_preds) %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~sample_type,scales = 'free')

full_join(epi_preds,endo_preds) %>% 
  ggplot(aes(x=dist_from_edge,y=observed)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~sample_type,scales = 'free')



endo_splines <- gdm::isplineExtract(gdm_endo) %>% as.data.frame() %>% mutate(sample_type='Endophyte')
epi_splines <- gdm::isplineExtract(gdm_epi) %>% as.data.frame() %>% mutate(sample_type='Epiphyte')

isplines <- full_join(endo_splines,epi_splines)
names(isplines) <- c("geographic_actual","dist_from_edge_actual",
                     "geographic_partial","dist_from_edge_partial","sample_type")

p1 <- isplines %>% 
  ggplot(aes(x=dist_from_edge_actual,y=dist_from_edge_partial,color=sample_type)) +
  geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2) +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        axis.text = element_text(face='bold',size=12)) +
  scale_color_manual(values=sampletypecolors) +
  labs(y="f(Distance from edge)",x="Distance from edge (m)",color="Sample type")


p2 <- isplines %>% 
  ggplot(aes(x=geographic_actual*100000,y=geographic_partial,color=sample_type)) +
  geom_smooth(se=FALSE,linewidth=2) +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        axis.text = element_text(face='bold',size=12)) +
  scale_color_manual(values=sampletypecolors) +
  labs(y="f(Geographic distance)",x="Geographic distance (m)",color="Sample type")

p1 + p2 + plot_layout(guides='collect')

ggsave("./Output/figs/Figure_4.png",dpi=300,height = 4,width = 10)
ggsave("./Output/figs/Figure_4.tiff",dpi=300,height = 4,width = 10)
