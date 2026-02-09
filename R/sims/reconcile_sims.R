library(tidyverse)

load(here::here("data/ovlp_sims1.rdata"))
load(here::here("data/ovlp_sims2.rdata"))
load(here::here("data/ovlp_sims3.rdata"))
load(here::here("data/ovlp_sims4.rdata"))


sim_out_diff_v2 <-
  bind_rows(
    sim_out_diff_v2.1,
    sim_out_diff_v2.2 %>%
      mutate(sim = sim + 250),
    sim_out_diff_v2.3 %>%
      mutate(sim = sim + 500),
    sim_out_diff_v2.4 %>%
      mutate(sim = sim + 750))


dist_comps <-
  bind_rows(
    dist_comps.1,
    dist_comps.2 %>%
      mutate(sim = sim + 250),
    dist_comps.3 %>%
      mutate(sim = sim + 500),
    dist_comps.4 %>%
      mutate(sim = sim + 750)
  )

cog_out_diff <-
  bind_rows(
    cog_out_diff.1,
    cog_out_diff.2,
    cog_out_diff.3,
    cog_out_diff.4
  )

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

save(sim_out_diff_v2,
     dist_comps,
     ovlp_ts_v2,
     file = here::here("data/ovlp_sims.rdata"))
