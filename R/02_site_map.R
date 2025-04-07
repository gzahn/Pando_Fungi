
# SETUP ####

# packages
if(!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
if(!requireNamespace("zahntools",quietly = TRUE)){devtools::install_github("gzahn/zahntools")}
library(tidyverse)
library(ggmap)
library(zahntools)

# Load google maps API key from .Renviron 
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private

# customization
mapstyle <- rjson::fromJSON(file = "./R/mapstyle.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

# data
path <- "./Data/Combined_Metadata_2025.csv"
dat <- read_csv(path) %>% 
  mutate(year = factor(collection_round))

latlon <- 
dat %>% 
  select(tree,lon,lat,year) %>% 
  filter(!is.na(lon) & !is.na(lat))

# BUILD MAP ####
area <- 
  ggmap::get_googlemap(center = c(lon = mean(latlon$lon), 
                                  lat = mean(latlon$lat)),
                       zoom = 15,
                       scale = 2,
                       style=mapstyle)
map <- 
ggmap::ggmap(area) +
  geom_point(data=latlon,aes(x=lon,y=lat,text=tree,color=year),size=2,alpha=1,shape=4) +
  scale_color_manual(values=c("red","black"))
map

plotly::ggplotly(map)

ggsave("./Output/site_map.png",dpi=300,height = 5,width = 5)
# file.copy("./Output/site_map.png","../gzahn.github.io/microbiome_bootcamp/media/site_map.png",overwrite = TRUE)
