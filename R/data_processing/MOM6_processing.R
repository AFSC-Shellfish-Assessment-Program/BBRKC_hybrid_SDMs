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


## Extract gridded data to reflect BT during tag release-------`
fns <- list.files("c:/users/seanh/documents/git/new_mom6/")
bt_all <- NULL
for (f in seq_along(fns)){
  print(fns[f])
  bt_init <- stars::read_ncdf(here::here("c:/users/seanh/documents/git/new_mom6/",fns[f]),
                              var = "tob", proxy = F)
  if (is.null(bt_all)) {
    bt_all <- bt_init
  } else {
    bt_all <- c(bt_all, bt_init)  # concatenate stars objects
  }

}

# bt_all
attr(bt_all,"dimensions")[[3]]$values <- seq.Date(from = as.Date("2005-01-01 12:00:00 UTC"),
                                                  by = "1 day", length.out = 7486)
# read MOM6 bottom temp data----
ocean_proj <- read_ncdf(here::here("data/mom6nep_hc202507_ocean_static_ak.nc"),
                        var = c("geolon","geolat"))

xc = matrix(ocean_proj[[1]], 342, 297)
yc = matrix(ocean_proj[[2]], 342, 297)

years <- 2005:2024
agg_temp_interannual_sum_aut <- NULL
big_grid_4326 <- big_grid %>% st_union() %>% st_transform(., "EPSG:4326") %>% st_make_valid()
for (i in years){
  for (m in c(10, 6)){
    message(paste0(i, m))

    # extract index for Oct 8-14 in each year
    if (m == 10){
      bt_ind_which <-
        which((year(as.POSIXct(attr(bt_all,"dimensions")[[3]]$values)) == i) &
                (month(as.POSIXct(attr(bt_all,"dimensions")[[3]]$values)) == m))[8:14]
    } else if (m == 6){
      bt_ind_which <-
        which((year(as.POSIXct(attr(bt_all,"dimensions")[[3]]$values)) == i) &
                (month(as.POSIXct(attr(bt_all,"dimensions")[[3]]$values)) == m))[1:7]
    }


    bt_slice <- bt_all[ , , , bt_ind_which] %>%
      st_as_stars(.,
                  curvilinear = list(ih = xc,
                                     jh = yc)) %>%
      .[big_grid_4326] %>%
      st_apply(., c("ih","jh"), mean)

    # convert to sf. Transform CRS to AEA coordinates
    sf_temp <-
      bt_slice %>%
      st_as_sf() %>%
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
