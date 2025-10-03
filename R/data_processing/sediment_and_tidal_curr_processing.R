library(stars)
library(tidyverse)
library(sf)
library(lubridate)
library(tsibble)
library(here)
library(raster)

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

agg_phi <- aggregate(sf_phi, by = big_grid, FUN = mean, na.rm = F) %>%
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

agg_tc <- aggregate(sf_tc, by = big_grid, FUN = mean, na.rm = FALSE) %>%
  st_difference(., ak_land) %>%
  dplyr::rename(tidal_curr = X3.pred)


# for (i in 1:nrow(agg_tc)){
#   tci <- agg_tc %>% slice(i) %>% pull(tidal_curr)
#   if (!is.na(tci)){
#     next
#   } else {
#     tc_buff <- agg_tc %>%
#       slice(i) %>%
#       st_buffer(dist = 10) %>%
#       dplyr::select(-tidal_curr) %>%
#       st_intersection(., agg_tc)
#
#     tc_buff_mean <- mean(tc_buff$tidal_curr, na.rm = T)
#     agg_tc[i,]$tidal_curr <- tc_buff_mean
#   }
# }
#
#
# ggplot( ) +
#   geom_sf(data = agg_tc, aes(fill = tidal_curr)) +
#   geom_sf(data = sf_s1) +
#   geom_sf(data = sf_s0)

save(agg_tc,
     file = here::here('data/agg_tc.rdata'))
