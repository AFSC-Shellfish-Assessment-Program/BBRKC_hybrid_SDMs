library(tidyverse)
library(raster)
library(sf)
library(stars)
library(Matrix)
library(mgcv)
library(sdmTMB)
library(RTMB)

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
  filter(!is.na(depth) & !is.na(temp) & !is.na(phi) & !is.na(tidal_curr)) %>%
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
                    year))

## deployment locations-----
sf_s0 = st_as_sf(male_crab_movement[,c('deploy_lon','deploy_lat','year')],
                 coords = c('deploy_lon','deploy_lat'),
                 crs = ak_crs)
s0 <- st_coordinates(sf_s0)

## pop-up locations-----
sf_s1 <- st_as_sf(male_crab_movement[,c('lon0','lat0','year')],
                  coords = c('lon0','lat0'),
                  crs = ak_crs)
s1 <- st_coordinates(sf_s1)

## deployment location grid cells----
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
## pop-up location grid cells----
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

## environmental data for fitting the model----
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

# load fitted models to get habitat preference formulas
load(here::here("data/fitted_movement_models.rdata"))

cs <- 25

# function to do LOOCV------
mm_loocv <- function(formula  = ~ 0 + tidal_curr + temp + depth,
                     data = out_interannual,
                     gridded_domain = grid2,
                     deployment_locs = grid0,
                     release_locs = grid1,
                     time = "year",
                     move_comps = c("diffusion", "taxis"),
                     apply_cov_scaling = FALSE,
                     d_scaling = "none",
                     cellsize = cs){

  gridded_domain_df <- gridded_domain %>%
    st_as_sf() %>%
    st_centroid() %>%
    dream::sfc_as_cols(., names = c("lon","lat"))

  pred_out <- NULL
  for (i in seq_along(deployment_locs)){
    message(i)
    grid0.loo <- deployment_locs[-i]
    grid1.loo <- release_locs[-i]

    if (i <= 13){
      tags_per_step = c(12, 13, 37)
      yr <- 2021
    } else if (i > 13 & i <= 26){
      tags_per_step = c(13, 12, 37)
      yr <- 2022
    } else if (i > 26){
      tags_per_step = c(13, 13, 36)
      yr <- 2023
    }

    mt_loo <-
      mm(formula,
         data,
         gridded_domain,
         deployment_locs = grid0.loo,
         release_locs = grid1.loo,
         time,
         tags_per_step,
         move_comps,
         cellsize,
         d_scaling,
         apply_cov_scaling)

    if (any(is.na(mt_loo$sd$par.fixed)) |
        any(is.na(mt_loo$sd$cov.fixed))){
      message(formula)
      stop()
    }

    if (all(c("diffusion","taxis") %in% move_comps)) {
      X_gk <- predict(mt_loo$preference_model,
                      newdata = out_interannual %>% filter(year == yr),
                      type = "lpmatrix")
      beta_k <- mt_loo$opt$par[names(mt_loo$opt$par) != "ln_D"]
      ln_D    <- mt_loo$opt$par[names(mt_loo$opt$par) == "ln_D"]
      h_s <- X_gk %*% beta_k
    } else if (identical(move_comps, "diffusion")) {
      ln_D    <- mt_loo$opt$par[names(mt_loo$opt$par) == "ln_D"]
      pref_g <- NULL
    }

    A_ss <- mt_loo$A
    At_zz = cbind( attr(A,"i"), attr(A,"j") ) + 1
    Mrate_gg <- M_dot(A_ss,
                      At_zz,
                      rate_par = ln_D,
                      pref_g = h_s,
                      move_comps,
                      d_scaling,
                      delta_d = cellsize)

    v_ig <- matrix(0, nrow=length(grid0), ncol=ncol(A))
    v_ig[cbind(seq_along(grid0), grid0)] = 1
    f_ig <- v_ig %*% t(M)

    v <- matrix(0, nrow = nrow(A), ncol = 1)
    v[grid0[i], 1] <- 1 # starting position for left out tag

    end_dist <- as.vector(t(Matrix::expm(Mrate_gg)) %*% v)
    # gridded_domain %>%
    #   st_as_sf() %>%
    #   mutate(end_dist) %>%
    #   ggplot() +
    #     geom_sf(aes(fill = end_dist)) +
    #     geom_sf(data = sf_s1 %>% slice(i))

    pred_int <- tibble(
      lon = weighted.mean(gridded_domain_df$lon, w = end_dist),
      lat = weighted.mean(gridded_domain_df$lat, w = end_dist),
      dropped_obs = i
    )
    assign("pred_out", rbind(pred_out, pred_int))
  }

  return(
    list(cv_preds = pred_out)
  )

}

m1_cv <-
  mm_loocv(formula  = m1.t$formula,
         data = out_interannual,
         gridded_domain = grid2,
         deployment_locs = grid0,
         release_locs = grid1,
         time = "year",
         move_comps = c("diffusion", "taxis"))

m2_cv <-
  mm_loocv(formula  = m2.t$formula,
           data = out_interannual,
           gridded_domain = grid2,
           deployment_locs = grid0,
           release_locs = grid1,
           time = "year",
           move_comps = c("diffusion", "taxis"))

m3_cv <-
  mm_loocv(formula  = m3.t$formula,
           data = out_interannual,
           gridded_domain = grid2,
           deployment_locs = grid0,
           release_locs = grid1,
           time = "year",
           move_comps = c("diffusion", "taxis"))

m4_cv <-
  mm_loocv(formula  = m4.t$formula,
           data = out_interannual,
           gridded_domain = grid2,
           deployment_locs = grid0,
           release_locs = grid1,
           time = "year",
           move_comps = c("diffusion", "taxis"))

m5_cv <-
  mm_loocv(formula  = m5.t$formula,
           data = out_interannual,
           gridded_domain = grid2,
           deployment_locs = grid0,
           release_locs = grid1,
           time = "year",
           move_comps = c("diffusion", "taxis"))

m6_cv <-
  mm_loocv(formula  = NULL,
           data = out_interannual,
           gridded_domain = grid2,
           deployment_locs = grid0,
           release_locs = grid1,
           time = "year",
           move_comps = "diffusion")

# function to calculate RMSE
# ~ 0 + te(depth, tidal_curr, k = 3) + s(temp, k = 4)
calc_error <- function(model_cv,
                       observed){
  spat_resid <-
    as.numeric(
      st_distance(
        model_cv$cv_preds %>%
          st_as_sf(., coords = c("lon", "lat"),
                   crs = ak_crs),
        observed,
        by_element = TRUE
      )
    )

  rmse <-
    data.frame(value = sqrt(mean(spat_resid^2)),
             metric = "RMSE")
  resid_df <- tibble(spat_resid)
  return(list(rmse = rmse,
              resid_df = resid_df))
}

movement_rmse <-
  bind_rows(
    calc_error(model_cv = m1_cv,
               observed = sf_s1)[[1]] %>%
      mutate(mod_num = 1),
    calc_error(model_cv = m2_cv,
               observed = sf_s1)[[1]] %>%
      mutate(mod_num = 2),
    calc_error(model_cv = m3_cv,
               observed = sf_s1)[[1]] %>%
      mutate(mod_num = 3),
    calc_error(model_cv = m4_cv,
               observed = sf_s1)[[1]] %>%
      mutate(mod_num = 4),
    calc_error(model_cv = m5_cv,
               observed = sf_s1)[[1]] %>%
      mutate(mod_num = 5),
    calc_error(model_cv = m6_cv,
               observed = sf_s1)[[1]] %>%
      mutate(mod_num = 6)
  )



movement_resids <-
  bind_rows(
    calc_error(model_cv = m1_cv,
               observed = sf_s1)[[2]] %>%
      mutate(mod_num = 1),
    calc_error(model_cv = m2_cv,
               observed = sf_s1)[[2]] %>%
      mutate(mod_num = 2),
    calc_error(model_cv = m3_cv,
               observed = sf_s1)[[2]] %>%
      mutate(mod_num = 3),
    calc_error(model_cv = m4_cv,
               observed = sf_s1)[[2]] %>%
      mutate(mod_num = 4),
    calc_error(model_cv = m5_cv,
               observed = sf_s1)[[2]] %>%
      mutate(mod_num = 5),
    calc_error(model_cv = m6_cv,
               observed = sf_s1)[[2]] %>%
      mutate(mod_num = 6)
  )
#
#
# ggplot(movement_resids) +
#   geom_histogram(aes(spat_resid)) +
#   facet_wrap(~mod_num, ncol = 1)

save(movement_rmse,
     movement_resids,
     file = here::here("data/movement_model_LOOCV_stats.rdata"))

