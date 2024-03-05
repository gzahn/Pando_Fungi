# Alpha Diversity

# SETUP ####

library(phyloseq)
library(vegan)
library(tidyverse)

ps <- readRDS("./Output/clean_phyloseq_object.RDS")
ps


# richness
# shannon div
# at different taxonomic levels
# difference between endo/epi?

estimate_richness(ps, measures = c('Observed', 'Shannon', 'Simpson'))
#does more than estimate richness
  #wants a phyloseq and some measures
  #observed is # of asvs
#problem:
  #sister chromatids could be counted as different asvs


ps_species <- 
  tax_glom(ps, "Species")
estimate_richness(ps_species, measures = c('Observed', 'Shannon', 'Simpson'))

# aGLOMerate as a tax rank

# View(ps_species@tax_table)
# 
# View(ps_species@sam_data)

# Add diversity estimates to metadata as new columns
  # Both glommed and nonglommed


#this not the best way, stupid col names
  # Do not do it like this. $ in col names causes errors

# ps@sam_data$obs_richness <- 
#   estimate_richness(ps, measures = c('Observed'))
# 
# ps@sam_data$shan_richness <- 
#   estimate_richness(ps, measures = c('Shannon'))
# 
# ps@sam_data$simp_richness <- 
#   estimate_richness(ps, measures = c('Simpson'))
# 
# ps@sam_data$obs_richness_glom <- 
#   estimate_richness(ps_species, measures = c('Observed'))
# 
# ps@sam_data$shan_richness_glom <- 
#   estimate_richness(ps_species, measures = c('Shannon'))
# 
# ps@sam_data$simp_richness_glom <- 
#   estimate_richness(ps_species, measures = c('Simpson'))
# 
# View(ps_species@sam_data)

#this way probably better
alpha_asv <- estimate_richness(ps, measures = c('Observed', 'Shannon', 'Simpson'))
ps_species <- tax_glom(ps, "Species")
alpha_spp <- estimate_richness(ps_species, measures = c('Observed', 'Shannon', 'Simpson'))

ps@sam_data$asv_richness <- alpha_asv$Observed
ps@sam_data$asv_shannon <- alpha_asv$Shannon
ps@sam_data$asv_simpson <- alpha_asv$Simpson

ps@sam_data$spp_richness <- alpha_spp$Observed
ps@sam_data$spp_shannon <- alpha_spp$Shannon
ps@sam_data$spp_simpson <- alpha_spp$Simpson


#this is not a dataframe, but still a phyloseq object
View(ps@sam_data)

#this way is lower level way to force it to work, 
sam <- ps@sam_data %>% as("data.frame")
# class(sam)
# 
# cor(sam$distance_from_edge, sam$sample_type)

#this how we look at difference between two things:
# lm(sam$distance_from_edge~sam$sample_type)

# Model ####
  #div ~ sample_type * dist_from_edge
  #def want interaction so can see dif between epi and endo

# Models
glm(data = sam,
    formula = asv_richness ~ sample_type * distance_from_edge)

glm(data = sam,
    formula = asv_shannon ~ sample_type * distance_from_edge)

glm(data = sam,
    formula = asv_simpson ~ sample_type * distance_from_edge)

glm(data = sam,
    formula = spp_richness ~ sample_type * distance_from_edge)

glm(data = sam,
    formula = spp_shannon ~ sample_type * distance_from_edge)

glm(data = sam,
    formula = spp_simpson ~ sample_type * distance_from_edge)




# MAP ####

# packages
if(!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
if(!requireNamespace("zahntools",quietly = TRUE)){devtools::install_github("gzahn/zahntools")}
library(tidyverse)
library(ggmap)
library(zahntools)


# Load google maps API key from .Renviron 
ggmap::register_google(key = "APIKEY") # Key kept private

# customization
mapstyle <- rjson::fromJSON(file = "./R/mapstyle.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

# data
# path <- "./Data/pando_sample_metadata_filled.csv"
dat <- sam

#wait, need endos and epis
View(sam)
sam_endo <- filter(sam, sample_type == "endophyte" )
sam_epi <- filter(sam, sample_type == "epiphyte" )

dat_endo <- sam_endo
dat_epi <- sam_epi


# ENDOPHYTE MAPS ####

latlon_endo <- 
  dat_endo %>% 
  select(tree,lon,lat, asv_richness:spp_simpson) %>% 
  filter(!is.na(lon) & !is.na(lat))

# BUILD MAP ####
area <- 
  ggmap::get_googlemap(center = c(lon = mean(latlon_endo$lon), 
                                  lat = mean(latlon_endo$lat)),
                       zoom = 15,
                       scale = 2,
                       style=mapstyle,
                       key = "APIKEY") # Keep key secret

# asv_richness
map_asv_rich <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_endo,aes(x=lon,y=lat, color= asv_richness),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend")+
  ggtitle("ASV ENDO Richness")

map_asv_rich

#asv_shannon
map_asv_shan <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_endo,aes(x=lon,y=lat, color= asv_shannon),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("ASV ENDO Shannon")
  
map_asv_shan

#asv_simp
map_asv_simp <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_endo,aes(x=lon,y=lat, color= asv_simpson),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("ASV ENDO Simpson")

map_asv_simp

#spp_richness
map_spp_rich <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_endo,aes(x=lon,y=lat, color= spp_richness),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("SPP ENDO Richness")
  
map_spp_rich

#SPP_Shannon
map_spp_shan <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_endo,aes(x=lon,y=lat, color= spp_shannon),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("SPP ENDO Shannon")

map_spp_shan

#SPP_Simpson
map_spp_simp <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_endo,aes(x=lon,y=lat, color= spp_simpson),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("SPP ENDO Simpson")

map_spp_simp



# plotly::ggplotly(map)
# getwd()
# ggsave("./site_map.png",dpi=300,height = 5,width = 5)

# List of variables
variables <- c("asv_richness", "asv_shannon", "asv_simpson", "spp_richness", "spp_shannon", "spp_simpson")

# Iterate through each variable
for (var in variables) {
  map <- ggmap::ggmap(area) +
    geom_point(data = latlon_endo, aes(x = lon, y = lat, color = .data[[var]]),
               size = 2, alpha = 0.7, shape = 10, stroke = 1) +
    scale_color_gradient(low = "yellow", high = "red", guide = "legend") +
    ggtitle(paste(toupper(substr(var, 1, 3)), substr(var, 5, nchar(var)), sep = " "))  # Capitalize and format the title
  
  print(map)  # Print the map
  
  # Save the map
  ggsave(paste0("alpha_maps/site_map_", var, "_endo.png"), plot = map, dpi = 300, height = 5, width = 5)
}


# EPIPHYTE MAPS ####

latlon_epi <- 
  dat_epi %>% 
  select(tree,lon,lat, asv_richness:spp_simpson) %>% 
  filter(!is.na(lon) & !is.na(lat))

# BUILD MAP ####
area <- 
  ggmap::get_googlemap(center = c(lon = mean(latlon_epi$lon), 
                                  lat = mean(latlon_epi$lat)),
                       zoom = 15,
                       scale = 2,
                       style=mapstyle,
                       key = "APIKEY") # Keep key secret

# asv_richness
map_asv_rich <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_epi,aes(x=lon,y=lat, color= asv_richness),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend")+
  ggtitle("ASV EPI Richness")

map_asv_rich

#asv_shannon
map_asv_shan <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_epi,aes(x=lon,y=lat, color= asv_shannon),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("ASV EPI Shannon")

map_asv_shan

#asv_simp
map_asv_simp <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_epi,aes(x=lon,y=lat, color= asv_simpson),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("ASV EPI Simpson")

map_asv_simp

#spp_richness
map_spp_rich <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_epi,aes(x=lon,y=lat, color= spp_richness),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("SPP EPI Richness")

map_spp_rich

#SPP_Shannon
map_spp_shan <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_epi,aes(x=lon,y=lat, color= spp_shannon),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("SPP EPI Shannon")

map_spp_shan

#SPP_Simpson
map_spp_simp <- 
  ggmap::ggmap(area) +
  geom_point(data=latlon_epi,aes(x=lon,y=lat, color= spp_simpson),
             size=2,alpha=0.7 ,shape=10, stroke = 1) +
  scale_color_gradient(low="yellow", high="red", guide = "legend") +
  ggtitle("SPP EPI Simpson")

map_spp_simp



# plotly::ggplotly(map)
# getwd()
# ggsave("./site_map.png",dpi=300,height = 5,width = 5)

# List of variables
variables <- c("asv_richness", "asv_shannon", "asv_simpson", "spp_richness", "spp_shannon", "spp_simpson")

# Saving the Maps

# Iterate through each variable
# for (var in variables) {
#   map <- ggmap::ggmap(area) +
#     geom_point(data = latlon_endo, aes(x = lon, y = lat, color = .data[[var]]),
#                size = 2, alpha = 0.7, shape = 10, stroke = 1) +
#     scale_color_gradient(low = "yellow", high = "red", guide = "legend") +
#     ggtitle(paste(toupper(substr(var, 1, 3)), substr(var, 5, nchar(var)), sep = " "))  # Capitalize and format the title
#   
#   print(map)  # Print the map
#   
#   # Save the map
#   ggsave(paste0("alpha_maps/site_map_", var, "_epi.png"), plot = map, dpi = 300, height = 5, width = 5)
# }

