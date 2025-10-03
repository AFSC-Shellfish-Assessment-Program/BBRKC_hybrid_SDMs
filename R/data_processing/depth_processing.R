library(stars)
library(tidyverse)
library(sf)
library(lubridate)
library(tsibble)
library(here)


ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"


# boundary and grid data for processing
load(here::here("data/spatial_layers.rdata"))

# depth----
sf_depth <-
  read_stars(here::here("data/efh_bathy_1km.tif")) %>%
  st_transform(ak_crs) %>%
  .[merged_boundaries] %>%
  st_as_sf() %>%
  dplyr::rename(depth = 1) %>%
  mutate(depth = ifelse(depth < 0, NA, depth))

## aggregate depth to grid----
agg_depth <- aggregate(sf_depth, by = big_grid, FUN = mean, na.rm = FALSE) %>%
  st_difference(., ak_land)

save(agg_depth,
     file = here::here('data/agg_depth.rdata'))
