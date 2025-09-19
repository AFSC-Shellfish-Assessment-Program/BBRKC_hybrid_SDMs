library(stars)
library(tidyverse)
library(here)
library(sf)
library(lubridate)
library(tsibble)

ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# Extend the BB management area to include those crab tagged just to the north of boundary
ebs <- st_read(here::here("data/ebs_grid.kml")) %>%
  st_zm() %>%
  st_transform(., ak_crs)

bb_mng_sf <-
  st_read(here::here("data/sap_layers.gdb"),
          layer = "BristolBaySurveyStrata", quiet = T)  %>%
  dplyr::rename(geometry = Shape) %>%
  st_transform(., crs= ak_crs) %>%
  summarise()

# create grid from bounding box
merged_boundaries <- bb_mng_sf %>%
  add_row(., ebs %>%
            dplyr::select(geometry)) %>%
  st_union()
big_grid <- st_make_grid(merged_boundaries,
                     cellsize = c(25,25)) %>%
  st_intersection(., merged_boundaries)

mid_grid <- st_make_grid(ebs,
                         cellsize = c(25,25)) %>%
  st_intersection(., ebs)

# Land polygon
ak_land <- rnaturalearthhires::states10 %>%
  filter(name == "Alaska") %>%
  st_union() %>%
  st_transform(., ak_crs) %>%
  st_intersection(., st_read(here::here("data/EBS_land.kml")) %>%
                    st_zm() %>%
                    st_transform(., ak_crs))

# other closure areas
other_layers <- st_layers(here::here("data/Closure areas"))$name
other_sa <- NULL
for (i in 1:length(other_layers)){
  # print(i)
  layer_init <- st_read(here::here(paste0("data/Closure areas")), quiet = T,
                        layer = other_layers[i])  %>%
    st_transform(., crs = 4326) %>%
    dplyr::summarise() %>%
    mutate(name = other_layers[i])

  assign('other_sa', rbind(other_sa, layer_init))
}

# management boundaries----
rkcsa_layers <- c("RKCSA_sub", "RKCSA")
rkcsa_sa <- NULL
for (i in rkcsa_layers){
  # print(i)
  layer_init <- st_read(here::here(paste0("data/management_boundaries/",i,".shp")), quiet = T)  %>%
    st_transform(., crs =4326) %>%
    dplyr::summarise() %>%
    mutate(name = "RKCSA")

  assign('rkcsa_sa', rbind(rkcsa_sa, layer_init))
}

manage_sf <-
  rkcsa_sa %>%
  dplyr::select(name, geometry) %>%
  bind_rows(.,
            st_read(here::here("data/sap_layers.gdb"),
                    layer = "BristolBaySurveyStrata", quiet = T)  %>%
              st_transform(., crs =4326) %>%
              dplyr::rename(geometry = Shape) %>%
              mutate(name = "BB survey") %>%
              dplyr::select(geometry, name))

save(big_grid,
     mid_grid,
     ak_land,
     manage_sf,
     bb_mng_sf,
     merged_boundaries,
     rkcsa_sa,
     ebs,
     other_sa,
     file = here::here("data/spatial_layers.rdata"))
