library(tidyverse)
library(raster)
library(patchwork)
library(sf)
library(stars)
library(sdmTMB)
library(INLA)
library(concaveman)
library(crabpack)
library(smoothr)

# load data------------
ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# plot limits
ylims = c(53.75, 62)
xlims = c(-175.5, -157.5)

bb_ll <- st_bbox(c(
  xmin = xlims[1], ymin = ylims[1],
  xmax = xlims[2], ymax = ylims[2]),
  crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  st_transform(., ak_crs)

## boundary and grid data---
load(here::here("data/spatial_layers.rdata"))

## load temperature data (october bottom temps from MOM6)----
load(here::here('data/agg_temp_interannual_sum_aut_mom6.rdata'))

## load bathymetry data----
load(here::here("data/agg_depth.rdata"))

## load sediment grain size data----
load(here::here("data/agg_phi.rdata"))

## load maxmimum tidal current data----
load(here::here("data/agg_tc.rdata"))

## load movement model grid for visualization
load(here::here("data/movement_model_particulars.rdata"))

## helper functions
source(here::here("R/helpers.R"))

# prediction grid--------
pred_grid <-
  agg_temp_interannual_sum_aut %>%
  filter(year %in% 2005:2024,
         month == 10) %>%
  st_join(., agg_depth %>%
            st_centroid(.)) %>%
  st_join(., agg_tc %>%
            st_centroid(.)) %>%
  filter(!is.na(depth) & !is.na(temp) & !is.na(tidal_curr)) %>%
  st_transform(., crs = ak_crs) %>%
  st_intersection(ebs)

# crabpack query----
# download and process NMFS EBS BTS data and calculate CPUE using crabpack
species <- "RKC"
region1 <- "EBS"
years <- 2005:2024
channel <- "API"

specimen_data_ebs_rkc <- crabpack::get_specimen_data(species = species,
                                                     region = region1, # Eastern Bering Sea
                                                     years = years,
                                                     channel = channel)

# get haul data for (gear) temperature and depth
haul_df_rkc <- specimen_data_ebs_rkc$haul %>%
  rename_all(str_to_lower) %>%
  dplyr::rename(latitude = mid_latitude,
                longitude = mid_longitude) %>%
  dplyr::select(bottom_depth,
                temp = gear_temperature,
                longitude,
                latitude,
                year,
                region,
                station_id)

# calculate CPUE
df_sf <-
  calc_cpue(crab_data = specimen_data_ebs_rkc,
            species = species,
            region = region1,
            year = years,
            crab_category = "mature_male") %>%
  dplyr::rename_all(str_to_lower) %>%
  dplyr::rename(lon = longitude,
                lat = latitude) %>%
  left_join(.,
            haul_df_rkc) %>%
  mutate(cpue = cpue/3.429) %>% # convert CPUE from crab nm^-2 to crab km^-2
  st_as_sf(., coords = c("lon","lat"), crs = 4326) %>%
  st_transform(., crs = ak_crs)

# mesh creation----
survey_buffed <-
  df_sf %>%
    filter(cpue > 0) %>% #
      concaveman::concaveman(., concavity = 3) %>% #Draw concave hull around positive catches
      st_buffer(., dist = 40) # extend buffer

ggplot() +
  geom_sf(data = survey_buffed) +
  geom_sf(data = ebs, fill = "transparent") +
  geom_sf(data = df_sf) +
  geom_sf(data = df_sf %>% filter(cpue > 0), color = "red")

# as data frame
df <- df_sf %>%
  st_intersection(., survey_buffed) %>%
  {. ->> df_sf_inter} %>%
  sfc_as_cols(., names = c("lon","lat")) %>%
  st_set_geometry(NULL) %>%
  mutate(yf = factor(year)) %>%
  filter(year >= 2005)

survey_buffed_sp <-
  survey_buffed %>%
  as_Spatial()

max.edge <- 30
bound.outer <- 100

mesh <- INLA::inla.mesh.2d(loc.domain = st_coordinates(df_sf_inter),
                           max.edge = c(1, 5) * max.edge,
                           offset = c(max.edge, bound.outer),
                           cutoff = max.edge)

## convert for sdmTMB----
ebs_spde <- sdmTMB::make_mesh(data = df,
                              xy_cols = c("lon", "lat"),
                              mesh = mesh)

plot(ebs_spde)
ebs_spde$mesh$n

# model fitting----------------------------------------------------------------------
tictoc::tic()
## tweedie
m1_sdm_tw <- sdmTMB::sdmTMB(
  formula = cpue ~ 0 + as.factor(year),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  mesh = ebs_spde,
  family = sdmTMB::tweedie(),
  data = df
)
tictoc::toc()
sanity(m1_sdm_tw)

## delta lognormal
tictoc::tic()
m1_sdm_dl1 <- sdmTMB::sdmTMB(
  formula = cpue ~ 0 + as.factor(year),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  mesh = ebs_spde,
  anisotropy = TRUE,
  family = sdmTMB::delta_lognormal(),
  data = df
)
tictoc::toc()
sanity(m1_sdm_dl1)

# cross-validation-------------------------------------------------------------------
# tweedie
m1_sdm_tw_cv <- sdmTMB::sdmTMB_cv(
  k_folds = 8,
  formula = cpue ~ 0 + as.factor(year),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  mesh = ebs_spde,
  family = sdmTMB::tweedie(),
  data = df
)

## delta lognormal
m1_sdm_dl1_cv <- sdmTMB::sdmTMB_cv(
  k_folds = 8,
  formula = cpue ~ 0 + as.factor(year),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  mesh = ebs_spde,
  anisotropy = TRUE,
  family = sdmTMB::delta_lognormal(),
  data = df
)

m1_sdm_dl1_cv$sum_loglik
m1_sdm_tw_cv$sum_loglik

# residuals---------------------------------------------------------------------------

## tweedie-------

### DHARMa-------
set.seed(67)
dharm_resids_tw <- simulate(m1_sdm_tw, nsim = 1000, type = "mle-mvn")
dharm_resids_tw2 <- dharma_residuals(dharm_resids_tw, m1_sdm_tw, return_DHARMa = TRUE)
DHARMa::testResiduals(dharm_resids_tw2)

### quantile residuals------
r1 <- residuals(m1_sdm_tw, type = "mle-mvn")
qqnorm(r1);abline(0, 1)
ks.test(r1, pnorm)

## delta lognormal------

### DHARMa------
set.seed(67)
dharm_resids <- simulate(m1_sdm_dl1, nsim = 1000, type = "mle-mvn")
dharm_resids2 <- dharma_residuals(dharm_resids, m1_sdm_dl1, return_DHARMa = TRUE)
DHARMa::testResiduals(dharm_resids2)

### quantile residuals-----
r2 <- residuals(m1_sdm_dl1,1, type = "mle-mvn")
qqnorm(r2);abline(0, 1)
ks.test(r2, pnorm)

r3 <- residuals(m1_sdm_dl1,2, type = "mle-mvn")
qqnorm(r3);abline(0, 1)
ks.test(r3, pnorm)

### MCMC residuals-----------
process <- FALSE
if (process){

  # some MCMC issues but residuals look good otherwise

  #### presence/absence-------
  samps1 <- sdmTMBextra::predict_mle_mcmc(m1_sdm_dl1,
                                          model = 1,
                                          mcmc_iter = 1000,
                                          mcmc_warmup = 500)
  mcmc_res1 <- residuals(m1_sdm_dl1, model = 1,
                         type = "mle-mcmc",
                         mcmc_samples = samps1)

  #### positive CPUE component---------
  samps2 <- sdmTMBextra::predict_mle_mcmc(m1_sdm_dl1,
                                          model = 2,
                                          mcmc_iter = 1000,
                                          mcmc_warmup = 500)
  mcmc_res2 <- residuals(m1_sdm_dl1, model = 2,
                         type = "mle-mcmc",
                         mcmc_samples = samps2)

  save(mcmc_res1,samps1, file = here::here("data/mcmc_res_pres_m1_dl.rdata"))
  save(mcmc_res2,samps2, file = here::here("data/mcmc_res_pres_m2_dl.rdata"))
} else {
  load(here::here("data/mcmc_res_pres_m1_dl.rdata"))
  load(here::here("data/mcmc_res_pres_m2_dl.rdata"))
}
qqnorm(mcmc_res1);abline(0,1)
ks.test(mcmc_res1, pnorm)
qqnorm(mcmc_res2);abline(0,1)
ks.test(mcmc_res2, pnorm)

### mesh figure-----
ak_land2 <- ak_land %>%
  st_intersection(, bb_ll)

png(filename = here::here("figs/supp_map/ebs_mesh.png"),
    width = 6, height = 6,
    units = "in", res = 300)
plot(ak_land2,
     axes = TRUE,
     col = "grey",
     xlim = c(-1100, -400),
     xlab = "Eastings (km)",
     ylab =  "Northings (km)",
     ylim = c(500, 1400))
plot(ebs_spde, add = T)
dev.off()

# prediction grid for spatiotemporal and movement models----
pred_grid_sfc <- pred_grid %>%
  filter(year == 2023, month == 10) %>%
  st_as_sfc()

# prediction grid figure----
png(filename = here::here("figs/supp_map/pred_grid.png"),
    width = 6, height = 6,
    units = "in", res = 300)
plot(ak_land2,
     axes = TRUE,
     col = "grey",
     xlim = c(-1100, -400),
     xlab = "Eastings (km)",
     ylab =  "Northings (km)",
     ylim = c(500, 1400))
plot(pred_grid_sfc, add = T, border = "purple", lwd = 2)
dev.off()

# expand prediction grid by year-----
nd <- pred_grid %>%
  st_centroid() %>%
  dplyr::select(year) %>%
  sfc_as_cols(., names = c("lon","lat")) %>%
  st_set_geometry(NULL) %>%
  filter(year != 2020) # no survey in 2020

# make predictions--------
pred <- sdmTMB:::predict.sdmTMB(object = m1_sdm_dl1,
                                newdata = nd,
                                type = "response") %>%
  st_as_sf(., coords = c("lon","lat"), crs = ak_crs)

# predictions as sf----
pred_df_sf <-
  pred_grid %>%
  dplyr::select(geometry) %>%
  distinct() %>%
  st_join(., pred) #%>%

# save outputs----
save(pred,
     df,
     pred_df_sf,
     m1_sdm_dl1,
     pred_grid,
     nd,
     file = here::here("data/st_model_predictions_dl.rdata"))

