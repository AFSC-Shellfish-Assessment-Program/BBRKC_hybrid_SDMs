library(tidyverse)
library(readxl)
library(lubridate)

# logbook data (unavailable online)-----
ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units= km +no_defs"

lbpre2024 <- read_csv(here::here("data/rkc.logbook.clean.csv")) %>%
  dplyr::select(-1) %>%
  st_as_sf(., coords = c('longitude','latitude'), crs = 4326) %>%
  st_transform(., crs = ak_crs)

files <- list.files(here::here("data","Clean DFL data"))
lb2024 <-
  read.csv(here::here("data","Clean DFL data",files[6])) %>%
  rename_all(str_to_lower) %>%
  rename_all(str_remove_all, pattern = "\\.") %>%
  filter(pot_count >= 5 & pot_count < 100) %>%
  rowwise() %>%
  mutate(longitude = mean(c(begin_lon,
                            end_lon)),
         latitude = mean(c(begin_lat,
                           end_lat)),
         date = as.Date(haul_date, format = "%m/%d/%Y"),
         year = year(date)) %>%
  dplyr::select(date, year, catch_pp = cpue, longitude, latitude) %>%
  st_as_sf(., coords = c('longitude','latitude'), crs = 4326) %>%
  st_transform(., crs = ak_crs)

lb <- bind_rows(lbpre2024,
                lb2024)

save(lb,
     file = here::here("data/dfl_data.rdata"))
