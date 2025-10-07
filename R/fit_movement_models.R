library(tidyverse)
library(raster)
library(sf)
library(stars)
library(Matrix)
library(mgcv)
library(sdmTMB)
library(RTMB)
library(splines2)

# Load data, source movement model----

## CRS----
ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units= km +no_defs"

## boundary and grid data----
load(here::here("data/spatial_layers.rdata"))

## load temperature layer (June and Otober bottom temps from MOM6)----
load(here::here('data/agg_temp_interannual_sum_aut_mom6.rdata'))

## load bathymetry layer----
load(here::here("data/agg_depth.rdata"))

## load sediment grain size layer----
load(here::here("data/agg_phi.rdata"))

## load tidal current maxima layer----
load(here::here("data/agg_tc.rdata"))

## source movement model----
source(here::here("R/mm.R"))

# helper functions
source(here::here("R/helpers.R"))

# Data processing----

## reduces size of spatial domain for model fitting----
ebs_abridged <- st_read(here::here("data/ebs_pretty_no_goa_move_lim.kml")) %>%
  st_zm() %>%
  st_transform(., crs = ak_crs)

## environmental covariates----
env_cov_full <- agg_temp_interannual_sum_aut %>%
  filter(year %in% 2005:2023,
         month == 10) %>%
  st_join(., agg_phi %>%
            st_centroid(.)) %>%
  st_join(., agg_depth %>%
            st_centroid(.)) %>%
  st_join(., agg_tc %>%
            st_centroid(.)) %>%
  filter(!is.na(depth) & !is.na(temp) & !is.na(tidal_curr) & !is.na(phi)) %>%
  st_transform(., crs = ak_crs) %>%
  st_intersection(., ebs_abridged)

env_cov <- env_cov_full %>% filter(year %in% 2021:2023,
                                   month == 10)

## import, process, and combine tagging data------
m2021 <-
  read_csv(here::here("data/2021_Oct_BBRKC_Male_FINAL.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = rel_lon,
         deploy_lat = rel_lat) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("deploy_lon","deploy_lat"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("deploy_lon","deploy_lat")) %>%
  st_set_geometry(NULL) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("lon0","lat0")) %>%
  st_set_geometry(NULL) %>%
  mutate(year = 2021)

m2022 <-
  read_csv(here::here("data/2022_Oct_BBRKC_Male_FINAL.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = rellongdd,
         deploy_lat = rellatdd) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("deploy_lon","deploy_lat"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("deploy_lon","deploy_lat")) %>%
  st_set_geometry(NULL) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("lon0","lat0")) %>%
  st_set_geometry(NULL) %>%
  mutate(year = 2022)

m2023 <-
  read_csv(here::here("data/2023_Oct_BBRKC_Male_FINAL.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = rellongdd,
         deploy_lat = releaselatdd) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("deploy_lon","deploy_lat"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("deploy_lon","deploy_lat")) %>%
  st_set_geometry(NULL) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("lon0","lat0")) %>%
  st_set_geometry(NULL) %>%
  mutate(year = 2023)

male_crab_movement <-
  bind_rows(
    m2021 %>%
      dplyr::select(tag, deploy_t, rel_t, deploy_days,
                    t0, lat0, lon0, deploy_lon, deploy_lat,
                    year) ,
    m2022 %>%
      dplyr::select(tag, deploy_t, rel_t, deploy_days,
                    t0, lat0, lon0, deploy_lon, deploy_lat,
                    year) ,
    m2023 %>%
      dplyr::select(tag, deploy_t, rel_t, deploy_days,
                    t0, lat0, lon0, deploy_lon, deploy_lat,
                    year)
  )

# deployment locations-----
sf_s0 = st_as_sf(male_crab_movement[,c('deploy_lon','deploy_lat','year')],
                 coords = c('deploy_lon','deploy_lat'),
                 crs = ak_crs)
s0 <- st_coordinates(sf_s0)

# pop-up locations-----
sf_s1 <- st_as_sf(male_crab_movement[,c('lon0','lat0','year')],
                  coords = c('lon0','lat0'),
                  crs = ak_crs)
s1 <- st_coordinates(sf_s1)

# deployment location grid cells----
grid0 <- c(unlist(st_intersects( sf_s0 %>%
                                   filter(year == 2021),
                                 env_cov %>%
                                   filter(year == 2021,
                                          month == 10))),
           unlist(st_intersects( sf_s0 %>%
                                   filter(year == 2022),
                                 env_cov %>%
                                   filter(year == 2022,
                                          month == 10))),
           unlist(st_intersects( sf_s0 %>%
                                   filter(year == 2023),
                                 env_cov %>%
                                   filter(year == 2023,
                                          month == 10))))

# pop-up location grid cells----
grid1 <- c(unlist(st_intersects( sf_s1 %>%
                                   filter(year == 2021),
                                 env_cov %>%
                                   filter(year == 2021,
                                          month == 10))),
           unlist(st_intersects( sf_s1 %>%
                                   filter(year == 2022),
                                 env_cov %>%
                                   filter(year == 2022,
                                          month == 10))),
           unlist(st_intersects( sf_s1 %>%
                                   filter(year == 2023),
                                 env_cov %>%
                                   filter(year == 2023,
                                          month == 10))))

# join spatial data to fit model
out_interannual <-
  env_cov %>%
  filter(month == 10) %>%
  group_by(year) %>%
  mutate(id = 1:n()) %>%
  st_set_geometry(NULL) %>%
  mutate(y = 0) %>%
  arrange(year, id) %>%
  ungroup() %>%
  mutate(id2 = 1:nrow(.),
         yf = factor(year))

## gridded spatial domain----
grid2 <-
  env_cov %>%
  filter(year == 2023) %>%
  st_as_sfc()

# cell size
cs <- 25

# Model fitting-----
m1.t <-
  mm(formula  = ~ 0 + tidal_curr + temp + depth,
     data = out_interannual,
     gridded_domain = grid2,
     deployment_locs = grid0,
     release_locs = grid1,
     time = "year",
     tags_per_step = c(13, 13, 37),
     move_comps = c("diffusion", "taxis"),
     cellsize = cs,
     d_scaling = "none")

m1.t$sd
m1.t$opt$message

m2.t <-
  mm(formula  = ~ 0 + s(tidal_curr, k = 4) + s(temp, k = 4) + s(depth, k = 4),
     data = out_interannual,
     gridded_domain = grid2,
     deployment_locs = grid0,
     release_locs = grid1,
     time = "year",
     tags_per_step = c(13, 13, 37),
     move_comps = c("diffusion", "taxis"),
     cellsize = cs,
     d_scaling = "none")

m2.t$sd
m2.t$opt$message

m3.t <-
  mm(formula  = ~ 0 + te(tidal_curr, temp, k = 3) + s(depth, k = 4),
     data = out_interannual,
     gridded_domain = grid2,
     deployment_locs = grid0,
     release_locs = grid1,
     time = "year",
     tags_per_step = c(13, 13, 37),
     move_comps = c("diffusion", "taxis"),
     cellsize = cs,
     d_scaling = "none")

m3.t$sd
m3.t$opt$message

m4.t <-
  mm(formula  = ~ 0 + te(depth, temp, k = 3) + s(tidal_curr, k = 4),
     data = out_interannual,
     gridded_domain = grid2,
     deployment_locs = grid0,
     release_locs = grid1,
     time = "year",
     tags_per_step = c(13, 13, 37),
     move_comps = c("diffusion", "taxis"),
     cellsize = cs,
     d_scaling = "none")

m4.t$sd
m4.t$opt$message

m5.t <-
  mm(formula  = ~ 0 + te(depth, tidal_curr, k = 3) + s(temp, k = 4),
     data = out_interannual,
     gridded_domain = grid2,
     deployment_locs = grid0,
     release_locs = grid1,
     time = "year",
     tags_per_step = c(13, 13, 37),
     move_comps = c("diffusion", "taxis"),
     cellsize = cs,
     d_scaling = "none")

m5.t$sd
m5.t$opt$message

# diffusion-only model
m6.t <-
  mm(formula  = NULL,
     data = out_interannual,
     gridded_domain = grid2,
     deployment_locs = grid0,
     release_locs = grid1,
     time = "year",
     tags_per_step = c(13, 13, 37),
     move_comps = "diffusion",
     cellsize = cs,
     d_scaling = "none")

m6.t$sd
m6.t$opt$message

save(m1.t,
     m2.t,
     m3.t,
     m4.t,
     m5.t,
     m6.t,
     file = here::here("data/fitted_movement_models.rdata"))

save(grid2,
     male_crab_movement,
     file = here::here("data/movement_model_particulars.rdata"))
