# Elevation Data

library(sf)
library(elevatr)


# Data frame with latitude and longitude
lat_lon_list <- data.frame(lat = sam$lat, lon = sam$lon)

# Convert to an sf object
lat_lon_sf <- st_as_sf(lat_lon_list, coords = c("lon", "lat"), crs = 4326)

# Get elevation data
sam$elev_data <- get_elev_point(lat_lon_sf)


# Un-uglify Columns
sam <- sam %>% 
  mutate(elevation = elev_data$elevation)

sam <- sam %>% 
  mutate(elev_geometry = elev_data$geometry)

sam <- sam %>% 
  mutate(elev_units = elev_data$elev_units)

sam <- sam %>% 
  select( -elev_data)

View(sam)
summary(sam$elevation)
