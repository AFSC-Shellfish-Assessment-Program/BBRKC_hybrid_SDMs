library(stars)
library(tidyverse)
library(sf)
library(lubridate)
library(tsibble)
library(here)

ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# boundary and grid data for processing
load(here::here("data/spatial_layers.rdata"))

# ebs area
ebs_area <- merged_boundaries %>%
  st_transform(., crs = 4326)

# load movement data------
m2021_s <-
  read_csv(here::here("data/2021_Oct_BBRKC_Male_FINAL.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = rel_lon,
         deploy_lat = rel_lat) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  mutate(year = 2021)

m2022_s <-
  read_csv(here::here("data/2022_Oct_BBRKC_Male_FINAL.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = rellongdd,
         deploy_lat = rellatdd) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  mutate(year = 2022)

m2023_s <-
  read_csv(here::here("data/2023_Oct_BBRKC_Male_FINAL.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = rellongdd,
         deploy_lat = releaselatdd) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  mutate(year = 2023)

# bind_rows(m2021_s,
#           m2022_s,
#           m2023_s) %>%
#   mutate(t0 = as.Date(t0, format = "%m/%d/%Y")) %>%
#   group_by(year) %>%
#   dplyr::summarise(r1 = range(t0)[1],
#                    r2 = range(t0)[2],
#                    avg_days = mean(deploy_days))
#
# bind_rows(m2021_s,
#           m2022_s,
#           m2023_s) %>%
#   mutate(deploy_t = as.Date(deploy_t, format = "%m/%d/%Y"),
#          dep_mon = deploy_t) %>%
#   group_by(year) %>%
#   dplyr::summarise(r1 = range(t0)[1],
#                    r2 = range(t0)[2],
#                    median = median(dep_mon))
#
# bind_rows(m2021_s,
#           m2022_s,
#           m2023_s) %>%
#   mutate(t0 = as.Date(t0, format = "%m/%d/%Y"),
#          deploy_t = as.Date(deploy_t, format = "%m/%d/%Y"),
#          rel_mon = yday(t0),
#          dep_mon = yday(deploy_t)) %>%
#   # group_by(year) %>%
#   dplyr::summarise(r_dep1 = range(dep_mon)[1],
#                    r_dep2 = range(dep_mon)[2],
#                    m_rel = median(rel_mon),
#                    m_dep = median(dep_mon),
#                    m = abs(m_rel - m_dep)/7)

## Extract gridded data to reflect BT during tag release-------`
fns <- list.files(here::here("data/new_mom6/"))
bt_all <- NULL
for (f in seq_along(fns)){
  print(fns[f])
  bt_init <- stars::read_ncdf(here::here("data/new_mom6/",fns[f]),
                              var = "tob", proxy = F)
  if (is.null(bt_all)) {
    bt_all <- bt_init
  } else {
    bt_all <- c(bt_all, bt_init)  # concatenate stars objects
  }

}
# bt_all
attr(bt_all,"dimensions")[[3]]$values <- seq.Date(from = as.Date("2005-01-01 12:00:00 UTC"),
                                                  by = "1 day", length.out = 7305)
# read MOM6 bottom temp data----
ocean_proj <- read_ncdf(here::here("data/mom6nep_hc202507_ocean_static_ak.nc"),
                        var = c("geolon","geolat"))

xc = matrix(ocean_proj[[1]], 342, 297)
yc = matrix(ocean_proj[[2]], 342, 297)

years <- 2005:2024
agg_temp_interannual_sum_aut <- NULL
for (i in years){
  for (m in c(10, 6)){
    message(paste0(i, m))

    # extract index for Oct 8-14 in each year
    bt_ind_which <-
      which((year(as.POSIXct(attr(bt_all,"dimensions")[[3]]$values)) == i) &
              (month(as.POSIXct(attr(bt_all,"dimensions")[[3]]$values)) == m))[8:14]

    bt_slice <- bt_all[ , , , bt_ind_which] %>%
      st_as_stars(.,
                  curvilinear = list(ih = xc,
                                     jh = yc)) %>%
      .[ebs_area] %>%
      st_apply(., c("ih","jh"), mean)

    # convert to sf and intersect with survey area. Transform CRS to AEA coordinates
    sf_temp <-
      bt_slice %>%
      st_as_sf() %>%
      st_intersection(., ebs_area) %>%
      dplyr::rename(temp = 1) %>%
      st_transform(ak_crs)

    agg_temp <- aggregate(sf_temp, by = big_grid, FUN = mean, na.rm = FALSE) %>%
      st_difference(., ak_land) %>%
      mutate(year = i,
             month = m)

    assign("agg_temp_interannual_sum_aut", rbind(agg_temp, agg_temp_interannual_sum_aut))
  }
}

save(agg_temp_interannual_sum_aut,
     file = here::here('data/agg_temp_interannual_sum_aut_mom6.rdata'))
