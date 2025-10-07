library(tidyverse)
library(sf)
library(stars)
library(sdmTMB)
library(mvtnorm)
library(furrr)

# Load data, source movement model----

## CRS----
ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

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

## fitted spatiotemporal model outputs-----
load(here::here("data/st_model_predictions.rdata"))

## fitted movement models----
load(here::here("data/fitted_movement_models.rdata"))

## movement model particulars----
load(here::here("data/movement_model_particulars.rdata"))

## source movement model----
source(here::here("R/mm.R"))

# helper functions
source(here::here("R/helpers.R"))

# logbook data-----
lb <- read_csv(here::here("data/rkc.logbook.clean.csv")) %>%
  dplyr::select(-1) %>%
  st_as_sf(., coords = c('longitude','latitude'), crs = 4326) %>%
  st_transform(., crs = ak_crs)

# load fishery CPUE data (aggregated)
load(here::here("data/agg_harvests_and_projections.rdata"))
lb_agg2 <- lb_agg %>% mutate(catch_pp = ifelse(is.na(catch_pp), 0, catch_pp))

# deployment locations-----
sf_s0 = st_as_sf(male_crab_movement[,c('deploy_lon','deploy_lat','year',
                                       'tag')],
                 coords = c('deploy_lon','deploy_lat'),
                 crs = ak_crs)
s0 <- st_coordinates(sf_s0)

# pop-up locations-----
sf_s1 <- st_as_sf(male_crab_movement[,c('lon0','lat0','year',
                                        'tag')],
                  coords = c('lon0','lat0'),
                  crs = ak_crs)
s1 <- st_coordinates(sf_s1)

#combined environmental covariates
env_cov_full <- agg_temp_interannual_sum_aut %>%
  filter(year %in% 2005:2023,
         month == 10) %>%
  st_join(., agg_phi %>%
            st_centroid(.)) %>%
  st_join(., agg_depth %>%
            st_centroid(.)) %>%
  st_join(., agg_tc %>%
            st_centroid(.)) %>%
  filter(!is.na(depth) & !is.na(temp) & !is.na(phi) & !is.na(tidal_curr)) %>%
  st_intersection(., ebs)

# deployment locations-----
sf_s0_locs <- male_crab_movement %>%
  dplyr::select(deploy_lon, deploy_lat, year, tag) %>%
  st_as_sf(., coords = c("deploy_lon","deploy_lat"),
           crs = ak_crs)

sf_s1_locs <- male_crab_movement %>%
  dplyr::select(lon0, lat0, year, tag) %>%
  st_as_sf(., coords = c("lon0","lat0"),
           crs = ak_crs)

# for predicting from the model
env_cov_full_df <-
  env_cov_full %>%
  st_set_geometry(NULL) %>%
  filter(month == 10) %>%
  mutate(y = 0)

# gridded domain
env_cov_full_sf <-
  env_cov_full %>%
       filter(year == 2005,
               month == 10) %>%
      st_as_sfc() %>%
      st_as_sf()

# gridded predictions from st model
nd_coords <- st_coordinates(nd %>%
                          filter(year == 2005) %>%
                          st_as_sf(., coords = c("lon","lat"),
                                        crs = ak_crs))

# the selected habitat preference model and the diffusion-only model
pref_mod <- m5.t
diff_mod <- m6.t

## estimated parameters from selected model
ests <- pref_mod$sd$par.fixed
est_diff <- diff_mod$sd$par.fixed

## fitted habitat preference model
f_mod1 <- pref_mod$model

# extract covariance matrix from fitted diffusion-taxis and diffusion-only model
cov_mat <- pref_mod$sd$cov.fixed
cov_diff <- diff_mod$sd$cov.fixed

# simulate from fitted SDM
sdm_sim <- predict(m1_sdm, newdata = nd,
                   nsim = nsim, type = "response")

# simulate from fitted movement models
par_draws_pref <- rmvnorm(nsim, mean = ests, sigma = cov_mat)
par_draws_diff <- rmvnorm(nsim, mean = est_diff, sigma = cov_diff)

# create adjacency matrix
st_rook <- function(m, ...) st_relate(m, m, pattern="F***1****", ... )
grid_A <- st_rook(env_cov_full %>%
                    filter(year == 2005,
                           month == 10) %>%
                    st_as_sfc(),
                  sparse=TRUE )
A_big <- as(grid_A,"sparseMatrix") |> as("TsparseMatrix")
At_zz = cbind( attr(A_big,"i"), attr(A_big,"j") ) + 1

all_years <- 2005:2023
survey_years <- c(2005:2019, 2021:2023)
sim_out_diff <- NULL
cog_out_diff <- NULL
for (i in all_years){

    if (i %in% survey_years){
      X_sz <- predict(pref_mod$preference_model, newdata  = env_cov_full_df %>% dplyr::filter(year == i), type = "lpmatrix")
      pred_year_sim_all <- sdm_sim[row.names(sdm_sim) == i,]

      for (j in 1:nsim){
        print(j)
        pred_year_sim <- pred_year_sim_all[,j] # simulated summer density vector

        # selected hab pref params and hab pref vector--
        ln_D_a <- par_draws_pref[j,][names(par_draws_pref[j,]) == "ln_D"]
        beta_k_a <- par_draws_pref[j,][names(par_draws_pref[j,]) != "ln_D"]
        h_s_a <- X_sz %*% cbind(beta_k_a)

        # diffusion only param--
        ln_D_b <- par_draws_diff[j,]

        # movement matricess---
        ## selected model sims
        Mrate_sim <- M_dot(A_ss = A_big,
                           At_zz,
                           rate_par = ln_D_a,
                           pref_g = h_s_a,
                           move_comps = c("diffusion","taxis"),
                           d_scaling = "none")
        M_sim <- Matrix::expm(Mrate_sim)

        ## diffusion only sims
        Mrate_sim_diff <- M_dot(A_ss = A_big,
                                At_zz,
                                rate_par = ln_D_b,
                                move_comps = "diffusion",
                                d_scaling = "none")
        M_sim_diff <- Matrix::expm(Mrate_sim_diff)

        # For overlap
        ovlp_int <-
          tibble(oct_proj = as.vector(t(pred_year_sim) %*% M_sim),
                 oct_diff = as.vector(t(pred_year_sim) %*% M_sim_diff),
                 p_proj  = oct_proj / sum(oct_proj),
                 p_diff  = oct_diff / sum(oct_diff),
                 lon = nd_coords[,1],
                 lat = nd_coords[,2],
                 year = i,
                 sim = j)

        # For COGs
        cog_int <-
          tibble(
            cog_lat = weighted.mean(nd_coords[, 2], ovlp_int$oct_proj),
            cog_lon = weighted.mean(nd_coords[, 1], ovlp_int$oct_proj),
            cog_lat_diff = weighted.mean(nd_coords[, 2], ovlp_int$oct_diff),
            cog_lon_diff = weighted.mean(nd_coords[, 1], ovlp_int$oct_diff),
            year = i,
            sim = j)

        assign("sim_out_diff", rbind(ovlp_int, sim_out_diff))
        assign("cog_out_diff", rbind(cog_int, cog_out_diff))

      }

    } else if (i == 2020) {
        ovlp_int <-
          tibble(oct_proj = NA,
                 oct_diff = NA,
                 p_proj = NA,
                 p_diff = NA,
                 lon = nd_coords[,1],
                 lat = nd_coords[,2],
                 year = i,
                 sim = NA)

        # For COGs
        cog_int <-
          tibble(
            cog_lat  = NA,
            cog_lon  = NA,
            cog_lat_diff  = NA,
            cog_lon_diff  = NA,
            year = i,
            sim = NA)

        assign("sim_out_diff", rbind(ovlp_int, sim_out_diff))
        assign("cog_out_diff", rbind(cog_int, cog_out_diff))
    }
  }

sim_out_diff_v2 <- NULL
for (yr in all_years){
    int2 <-
      env_cov_full_sf %>%
      st_join(.,
              sim_out_diff %>%
                filter(year == yr) %>%
                st_as_sf(., coords = c("lon","lat"), crs= ak_crs)) %>%
      st_join(.,
              lb_agg2  %>%
                filter(year == yr) %>%
                st_centroid(.) %>%
                mutate(p_catch = catch_pp/sum(catch_pp)))  %>%
      st_set_geometry(NULL) %>%
      dplyr::rename(year = year.x) %>%
      dplyr::select(-year.y)

  assign("sim_out_diff_v2", rbind(int2, sim_out_diff_v2))
}


# we can calculate COGs in 2021/22 because the survey happened, but overlap cannot be calculated.
ovlp_ts_v2 <-
  sim_out_diff_v2 %>%
    group_by(year, sim) %>%
    dplyr::summarise(bhat_proj = sum(sqrt(p_catch * p_proj)),
                     bhat_diff = sum(sqrt(p_catch * p_diff))) %>%
    left_join(., cog_out_diff) %>%
    {. ->> dist_comps} %>%
    group_by(year) %>%
    dplyr::summarise(mean_proj = mean(bhat_proj, na.rm = T),
                     upr95_proj = quantile(bhat_proj, 0.975, na.rm = T),
                     lwr95_proj = quantile(bhat_proj, 0.025, na.rm = T),

                     mean_diff = mean(bhat_diff, na.rm = T),
                     upr95_diff = quantile(bhat_diff, 0.975, na.rm = T),
                     lwr95_diff = quantile(bhat_diff, 0.025, na.rm = T),

                     mean_cog_lat = mean(cog_lat, na.rm = T),
                     upr95_cog_lat = quantile(cog_lat, 0.975, na.rm = T),
                     lwr95_cog_lat = quantile(cog_lat, 0.025, na.rm = T),

                     mean_cog_lon = mean(cog_lon, na.rm = T),
                     upr95_cog_lon = quantile(cog_lon, 0.975, na.rm = T),
                     lwr95_cog_lon = quantile(cog_lon, 0.025, na.rm = T),


                     mean_cog_lat_diff = mean(cog_lat_diff, na.rm = T),
                     upr95_cog_lat_diff = quantile(cog_lat_diff, 0.975, na.rm = T),
                     lwr95_cog_lat_diff = quantile(cog_lat_diff, 0.025, na.rm = T),

                     mean_cog_lon_diff = mean(cog_lon_diff, na.rm = T),
                     upr95_cog_lon_diff = quantile(cog_lon_diff, 0.975, na.rm = T),
                     lwr95_cog_lon_diff = quantile(cog_lon_diff, 0.025, na.rm = T))
