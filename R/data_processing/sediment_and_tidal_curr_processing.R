library(stars)
library(tidyverse)
library(sf)
library(lubridate)
library(tsibble)
library(here)

ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# boundary and grid data for processing
load(here::here("data/spatial_layers.rdata"))

# sediment grain size-----
# not used in analysis
sf_phi <- raster(here::here("data/EBS_phi_1km.gri")) %>%
  st_as_stars() %>%
  st_transform(., crs = ak_crs) %>%
  .[merged_boundaries] %>%
  st_as_sf()

agg_phi <- aggregate(sf_phi, by = big_grid, FUN = mean, na.rm = T) %>%
  st_difference(., ak_land) %>%
  dplyr::rename(phi = X3.pred)

save(agg_phi,
     file = here::here('data/agg_phi.rdata'))

# tidal current maxima----
sf_tc <- raster(here::here("data/Tmax.gri")) %>%
  st_as_stars() %>%
  st_transform(., crs = ak_crs) %>%
  .[merged_boundaries] %>%
  st_as_sf()

agg_tc <- aggregate(sf_tc, by = big_grid, FUN = mean, na.rm = T) %>%
  st_difference(., ak_land) %>%
  dplyr::rename(tidal_curr = X3.pred)

save(agg_tc,
     file = here::here('data/agg_tc.rdata'))
