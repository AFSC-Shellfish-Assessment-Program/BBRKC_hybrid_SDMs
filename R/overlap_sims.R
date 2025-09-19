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

nsim <- 1000
DeltaT <- 19

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

# create adjacency matrix
st_rook <- function(m, ...) st_relate(m, m, pattern="F***1****", ... )
grid_A <- st_rook(env_cov_full %>%
                    filter(year == 2005,
                           month == 10) %>%
                    st_as_sfc(),
                  sparse=TRUE )
A_big <- as(grid_A,"sparseMatrix") |> as("TsparseMatrix")

# runs in parallel for speed
plan(multisession, workers = parallel::detectCores() - 4)

# sim constants
set.seed(20250418)
survey_years <- c(2005:2019, 2021:2023)
tag_years <- 2021:2023
fishing_years <- c(2005:2020, 2023)
all_years <- 2005:2023

# extract environmental data
out <- env_cov_full %>%
  st_set_geometry(NULL) %>%
  mutate(y = 0)

# load fishery CPUE data (aggregated)
load(here::here("data/agg_harvests_and_projections.rdata"))
lb_agg2 <- lb_agg %>% mutate(catch_pp = ifelse(is.na(catch_pp), 0, catch_pp))
nd_sf <- nd %>% st_as_sf(., coords = c("lon","lat"), crs = ak_crs)

# do basis expansion by year and add to list
lpmats <- list()
for (y in 1:length(survey_years)){
  int <- predict(pref_mod$preference_model,
                 newdata  = out %>% dplyr::filter(year == !!survey_years[y]),
                 type = "lpmatrix")
  lpmats[[y]] <- int
}
names(lpmats) <- as.character(survey_years)

# simulate from fitted SDM
sdm_sim <- predict(m1_sdm, newdata = nd,
                   nsim = nsim, type = "response")

# used to extract predictions from sdmTMB simulations later on
sdm_rows_by_year <- split(
  seq_len(nrow(nd)),
  nd$year             # year is grouping factor
)

# Tag‑endpoint helper (these two functions are specific to calculating predicted distance to observed pop-up locations, not projections)
# takes tag specific movement probabilities (prob_mat) and finds the cell
# with the highest likelihood of movement. returns the cell lon/lat, tag,
# and whether likelihood is from diffusion-taxis model or diffusion-only model
predict_tag_endpoint <- function(deployment,
                                 td_mat = M_sim,
                                 do_mat = M_sim_diff,
                                 proj_grid_template = env_cov_full_sf,
                                 sim = draw_idx,
                                 year = yr,
                                 adj = A_big) {

  grid0 <- unlist(st_intersects(deployment, proj_grid_template))
  v_ig <- matrix(0, nrow=length(grid0), ncol=ncol(adj))
  v_ig[cbind(seq_along(grid0), grid0)] = 1

  f_ig_sim  <- v_ig %*% t(td_mat)
  f_ig_diff <- v_ig %*% t(do_mat)

  # extract movement probabilities for individual tags
  prob <- NULL
  for (i in 1:nrow(deployment)){

    tag <- deployment %>% slice(i) %>% pull(tag)

    t_ind_td <- which.max(f_ig_sim[i,])
    max_prob_td <- f_ig_sim[i,t_ind_td]
    loc_td <- env_cov_full_sf %>% slice(t_ind_td) %>% st_centroid() %>% st_coordinates()

    t_ind_diff <- which.max(f_ig_diff[i,])
    max_prob_diff <- f_ig_diff[i,t_ind_diff]
    loc_diff <- env_cov_full_sf %>% slice(t_ind_diff) %>% st_centroid() %>% st_coordinates()

    assign("prob", rbind(prob, tibble(
      tag,
      grid_index_td = t_ind_td,
      max_prob_td,
      grid_index_diff = t_ind_diff,
      max_prob_diff,
      lon_td = loc_td[,"X"],
      lat_td = loc_td[,"Y"],
      lon_diff = loc_diff[,"X"],
      lat_diff = loc_diff[,"Y"],
      sim,
      year
        )
      )
    )
  }
  return(prob)
}

par_draws_best <- rmvnorm(nsim, mean = ests, sigma = cov_mat)
par_draws_diff <- rmvnorm(nsim, mean = est_diff, sigma = cov_diff)

At_zz = cbind( attr(A_big,"i"), attr(A_big,"j") ) + 1
colsumA_g <- colSums(A_big)

# calculate overlap for simulated projection from diffusion+taxis and diffusion-only model
run_single_sim <- function(draw_idx) {

    beta_tau <- par_draws_best[draw_idx, ] #get parameter draws
    diff_par <- par_draws_diff[draw_idx, ]

    ovlp_lst <- list()
    dist_lst <- list()

    for (yr in all_years) {

      if (yr %in% survey_years){
        Xgk    <- lpmats[[as.character(yr)]]
        nbeta  <- ncol(Xgk)
        beta_k <- beta_tau[names(beta_tau) != "ln_D"] # habitat preference parameters
        tau    <- beta_tau[names(beta_tau) == "ln_D"] # diffusion parameter from best fitting model
        h_s <- Xgk %*% cbind(beta_k)
        Mrate_sim <- make_M(A_ss = A_big,
                            At_zz,
                            rate_par = tau,
                            pref_g = h_s,
                            move_comps = c("diffusion","taxis"),
                            d_scaling = "scale2",
                            DeltaD = 25,
                            colsumA_g)
        M_sim <- Matrix::expm(DeltaT * Mrate_sim)

        Mrate_sim_diff <- make_M(A_ss = A_big,
                                  At_zz,
                                  rate_par = tau,
                                  pref_g = h_s,
                                  move_comps = "diffusion",
                                  d_scaling = "scale2",
                                  DeltaD = 25,
                                  colsumA_g)
        M_sim_diff <- Matrix::expm(DeltaT * Mrate_sim_diff)
        idx_rows  <- sdm_rows_by_year[[as.character(yr)]]
        sdm_slice <- sdm_sim[idx_rows, draw_idx] # pull out year-specific June predictions from sdmTMB simulations

        oct_proj <- as.vector(t(sdm_slice) %*% M_sim) # project using diffusion-taxis...
        oct_diff <- as.vector(t(sdm_slice) %*% M_sim_diff)# or diffusion-only model

        p_proj  <- oct_proj / sum(oct_proj)
        p_diff  <- oct_diff / sum(oct_diff)

        coords  <- st_coordinates(nd_sf %>% filter(year == yr)) # year-specific grid locations

        cog_lat <- weighted.mean(coords[, 2], oct_proj)
        cog_lon <- weighted.mean(coords[, 1], oct_proj)

        cog_lat_diff = weighted.mean(coords[, 2], oct_diff)
        cog_lon_diff = weighted.mean(coords[, 1], oct_diff)
      }

      if (yr %in% fishing_years){
        catch_rows <- lb_agg2 %>% filter(year == yr) # year-specific aggregated CPUE
        catches <- catch_rows$catch_pp
        p_catch <- catches / sum(catches) # proportions for calculating overlap
      }

      if ((yr %in% survey_years & yr %in% fishing_years)){
        bhat_proj    <- sum(sqrt(p_catch * p_proj))
        bhat_diff    <- sum(sqrt(p_catch * p_diff))
      } else if ((yr %in% survey_years & !yr %in% fishing_years)){
        bhat_proj    <- NA
        bhat_diff    <- NA
      } else if (!yr %in% survey_years & yr %in% fishing_years){
        bhat_proj    <- NA
        bhat_diff    <- NA
        cog_lat <- NA
        cog_lon <- NA
        cog_lat_diff <- NA
        cog_lon_diff <- NA
      }

      ovlp_lst[[length(ovlp_lst) + 1]] <- tibble(
        year = yr,
        sim = draw_idx,
        bhat_proj,
        bhat_diff,
        cog_lon_diff,
        cog_lat_diff,
        cog_lon,
        cog_lat
      )

      # tag‑validation runs regardless of catches in years
      # where tags were deployed
      if (yr %in% tag_years) {
        sf_s0_init <- sf_s0 %>% filter(year == yr)
        message("Year ", yr, ": n_tags = ", nrow(sf_s0_init))
        dist_tbl <- predict_tag_endpoint(deployment = sf_s0_init,
                                         td_mat = M_sim,
                                         do_mat = M_sim_diff,
                                         proj_grid_template = env_cov_full_sf,
                                         sim = draw_idx,
                                         year = yr,
                                         adj = A_big)
        dist_lst[[length(dist_lst) + 1]] <- dist_tbl
      }
    }
    # ---- parameter draws -----------------------------------------
    list(ovlp = bind_rows(ovlp_lst),
         dist = bind_rows(dist_lst))
}

# ---- parallel simulation -------------------------------------
res_lst <- future_map(1:nsim, run_single_sim,
                      .options = furrr_options(seed = TRUE))

sim_out_diff_v2 <- bind_rows(map(res_lst, "ovlp"))
dist_sim_out_v2 <- bind_rows(map(res_lst, "dist")) %>%
  left_join(.,
            male_crab_movement %>%
              dplyr::select(lon0, lat0, tag, year),
            by = c("year","tag"))

ovlp_ts_v2 <-
  sim_out_diff_v2 %>% # we can calculate COGs in 2021/22 because the survey happened, but overlap cannot be calculated.
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
                   lwr95_cog_lon_diff = quantile(cog_lon_diff, 0.025, na.rm = T)) %>%
  add_row(year = c(2020, 2021, 2022))


dist_sim_out_v2$td_distance <-
  as.numeric(
    st_distance(
      dist_sim_out_v2 %>%
        st_as_sf(., coords = c("lon_td", "lat_td"),
                 crs = ak_crs),

      dist_sim_out_v2  %>%
        st_as_sf(., coords = c("lon0", "lat0"),
                 crs = ak_crs),
      by_element = TRUE
    )
  )


dist_sim_out_v2$diff_distance <-
  as.numeric(
    st_distance(
      dist_sim_out_v2 %>%
        st_as_sf(., coords = c("lon_diff", "lat_diff"),
                 crs = ak_crs),

      dist_sim_out_v2  %>%
        st_as_sf(., coords = c("lon0", "lat0"),
                 crs = ak_crs),
      by_element = TRUE
    )
  )

# summarise predicted tag pop up locations by diffusion-taxis and diffusion-only models
dist_ts <-
dist_sim_out_v2 %>%
  dplyr::select(td_distance, diff_distance, year, tag) %>%
  gather(., model, distance, -year, -tag) %>%
  mutate(model = ifelse(model == "td_distance", "diffusion + taxis", "diffusion only")) %>%
  group_by(year, model, tag) %>%
  summarise(
    mean_dist = mean(distance),
    lwr95 = quantile(distance, 0.025),
    upr95 = quantile(distance, 0.975))

