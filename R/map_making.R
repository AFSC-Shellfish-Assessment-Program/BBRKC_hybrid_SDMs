library(tidyverse)
library(raster)
library(sf)
library(pals)
library(gganimate)
library(stars)
library(ggspatial)
library(ggnewscale)
library(patchwork)
library(smoothr)
library(ggtext)

# helper functions
source(here::here("R/helpers.R"))

# data processing----
# Extend the BB management area to include those crab tagged just to the north of boundary
load(here::here("data/spatial_layers.rdata"))

ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# create grid from bounding box
polygon_sf <-
  ebs %>%
  st_union(.,bb_mng_sf) %>%
  summarise()

# plot limits
ylims = c(54.6, 58.95)
xlims = c(-167.8, -157.5)

bb_ll <- st_bbox(c(
  xmin = xlims[1], ymin = ylims[1],
  xmax = xlims[2], ymax = ylims[2]),
    crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  st_transform(., ak_crs) %>%
  st_bbox()

# load temperature data---
load(here::here('data/agg_temp_interannual_sum_aut_mom6.rdata'))

# load bathymetry data----
load(here::here("data/agg_depth.rdata"))

# load maxmimum tidal current data----
load(here::here("data/agg_tc.rdata"))

env_cov_full <-
  agg_temp_interannual_sum_aut %>%
  filter(year %in% 2005:2024) %>%
  st_join(., agg_depth %>%
            st_centroid(.)) %>%
  st_join(., agg_tc %>%
            st_centroid(.)) %>%
  filter(!is.na(depth) & !is.na(temp) & !is.na(tidal_curr)) %>%
  st_transform(., crs = ak_crs) %>%
  st_intersection(ebs)

# process bathymetry data----
sf_depth <-
  read_stars(here::here("data/efh_bathy_1km.tif")) %>%
  st_warp(., crs = ak_crs) %>%
  {. ->> all_depth} %>%
  .[polygon_sf] %>%
  st_as_sf() %>%
  dplyr::rename(depth = 1) %>%
  mutate(depth = ifelse(depth < 0, NA, depth))

# 100 m and 50 m contours
dc100 <- all_depth %>%
  st_contour(x = .,
             breaks = 100,
             na.rm = T,
             contour_lines = T) %>%
  smoothr::drop_crumbs(., threshold = 25) %>%
  st_intersection(., sf_depth) %>%
  st_transform(., 4326)

dc75 <- all_depth %>%
  st_contour(x = .,
             breaks = 75,
             na.rm = T,
             contour_lines = T) %>%
  smoothr::drop_crumbs(., threshold = 25) %>%
  st_intersection(., sf_depth) %>%
  st_transform(., 4326)

dc50 <- all_depth %>%
  st_contour(x = .,
             breaks = 50,
             na.rm = T,
             contour_lines = T) %>%
  smoothr::drop_crumbs(., threshold = 25) %>%
  st_intersection(., sf_depth) %>%
  st_transform(., 4326)
rm(all_depth)

## crab data-----
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

m2024 <-
  read_csv(here::here("data/2024_Oct_BBRKC_Male_Final.csv")) %>%
  rename_all(str_to_lower) %>%
  rename_all(function(x){str_replace(x, "\\.", "_")}) %>%
  rename(deploy_lon = release_lon,
         deploy_lat = release_lat) %>%
  mutate(tag = factor(tag)) %>%
  st_as_sf(., coords = c("deploy_lon","deploy_lat"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("deploy_lon","deploy_lat")) %>%
  st_set_geometry(NULL) %>%
  st_as_sf(., coords = c("lon0","lat0"), crs = 4326) %>%
  st_transform(., crs = ak_crs) %>%
  sfc_as_cols(., names = c("lon0","lat0")) %>%
  st_set_geometry(NULL) %>%
  mutate(year = 2024) %>%
  filter(tag != "263793")

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
                    year),
    m2024 %>%
      dplyr::select(tag, deploy_t, rel_t, deploy_days,
                    t0, lat0, lon0, deploy_lon, deploy_lat,
                    year)
  )
## survey locations-----
surv_locs <-
  read_csv(file = here::here("data/RKC_matmale_matfem_2005-2024.csv")) %>%
  dplyr::rename_all(str_to_lower) %>%
  dplyr::rename(year = survey_year,
                lon = longitude,
                lat = latitude) %>%
  filter(year %in% 2023,
         mat_sex == "Mature Male") %>%
  st_as_sf(., coords = c("lon", "lat"),
           crs = 4326) %>%
  st_intersection(., polygon_sf %>%
                    st_transform(., 4326)) %>%
  st_transform(., ak_crs) %>%
  sfc_as_cols(., names = c('lon', 'lat')) %>%
  st_set_geometry(NULL)

# NBBTCA
nbbtca <-
  st_read(here::here("data/162W_boundary_line.kml")) %>%
    st_zm() #%>%

# map work----
n <- 20
depth_pal <- pals::kovesi.linear_blue_95_50_c20(n)
depth_pal <- rev(pals::ocean.ice(n))

depth_plt <-
  ggplot(agg_depth) +
  geom_sf(aes(fill = depth, color = depth)) +
  scale_fill_gradientn(colors = depth_pal) +
  scale_color_gradientn(colors = depth_pal) +
  geom_sf(data = dc100 %>% summarise(), color = "grey50", alpha = 1, linewidth = 0.35) +
  geom_sf(data = dc75 %>% summarise(), color = "grey50", alpha = 1, linewidth = 0.35) +
  geom_sf(data = dc50 %>% summarise(), color = "grey50", alpha = 1, linewidth = 0.35) +
  guides(fill = guide_colorbar(
    barheight = unit(0.25, "in"),
    barwidth  = unit(3, "in"),
    direction     = "horizontal",
    title.position = "left",
    title.vjust   = 0.8,
    title = "Depth (m)"),
    color = "none") +
  ggnewscale::new_scale_color() +
  geom_sf(data = manage_sf,
          aes(color = name),
          fill = "transparent",
          linewidth = 0.5,
          show.legend = F) +
  scale_color_manual(values = c("black","purple"))+
  ggnewscale::new_scale_color() +
  geom_sf(data = nbbtca, color = "purple",
          linewidth = 0.5,
          linetype= 2,
          show.legend = F,
          fill = "transparent") +
  # land---
  geom_sf(data = ak_land) +
  # survey locations---
  geom_point(data = surv_locs,
             aes(y = lat, x = lon),
             shape = "+",
             size = 2,
             alpha = 0.75) +
  # tagging data---
  ggnewscale::new_scale_color() +
  geom_point(data = male_crab_movement,
             aes(y = deploy_lat,
                 x = deploy_lon),
             size = 0.5,
             color = "grey10",
             alpha = 0.75,
             show.legend = F) +
  geom_segment(data = male_crab_movement,
               aes(y = deploy_lat, x = deploy_lon,
                   yend = lat0, xend = lon0),
               arrow = arrow(length = unit(0.05, "inches")), show.legend = F,
               color = "grey10",
               alpha = 0.75,
               linewidth = 0.25) +
  labs(y = "Latitude",
       x = "Longitude") +
  theme_bw() +
    coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
             ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
             expand = FALSE) +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(angle = 0))

ggsave(depth_plt,
       filename = here::here("figs/fig1_depth_plt.svg"),
       dpi = 300,
       height = 6,
       width = 6)

## bottom temp----
bwr2 <- pals::kovesi.diverging_bwr_40_95_c42(100)
bt_plt <-
agg_temp_interannual_sum_aut %>%
  mutate(month_name = ifelse(month == 6, "June", "October")) %>%
  mutate(year_f = factor(year)) %>%
  filter(year %in% 2021:2024) %>%
    ggplot() +
    geom_sf(aes(fill = temp, color = temp)) +

    # land---
    scale_fill_gradientn(colors=bwr2,
                         breaks = c(2.5, 5, 7.5, 10)) +
    scale_color_gradientn(colors=bwr2,
                       breaks = c(2.5, 5, 7.5, 10)) +

    guides(fill = guide_colorbar(
      barheight = unit(0.25, "in"),
      barwidth  = unit(3, "in"),
      direction     = "horizontal",
      title.position = "left",
      title.vjust   = 0.8,
      title = "Bottom temp. (°C)"
    ),
    color = "none") +

    # management boundaries---
  ggnewscale::new_scale_color() +
  geom_sf(data = manage_sf,
          aes(color = name),
          fill = "transparent",
          linewidth = 0.5,
          show.legend = F) +
  scale_color_manual(values = c("black","purple"))+
  ggnewscale::new_scale_color() +
  geom_sf(data = nbbtca, color = "purple",
          linewidth = 0.5,
          linetype= "11",
          show.legend = F,
          fill = "transparent") +

    # tagging data---
  geom_point(data = male_crab_movement %>%
               mutate(month = 6),
             aes(y = deploy_lat,
                 x = deploy_lon),
             size = 0.5,
             color = "grey10",
             alpha = 0.75,
             show.legend = F) +
  geom_point(data = male_crab_movement %>%
               mutate(month = 10),
             aes(y = lat0,
                 x = lon0),
             size = 0.5,
             color = "grey10",
             alpha = 0.75,
             show.legend = F) +
  geom_segment(data = male_crab_movement,
               aes(y = deploy_lat, x = deploy_lon,
                   yend = lat0, xend = lon0),
               color = "grey10",
               alpha = 0.75,
               linewidth = 0.25) +
    # land---
    geom_sf(data = ak_land) +
    # etc---
    facet_wrap(year~month, ncol = 2) +
    labs(y = "Latitude",
         x = "Longitude") +
    theme_bw() +
    coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
             ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
             expand = FALSE) +
    theme(strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA,
                                      size = 0.2),
          legend.key = element_blank(),
          axis.title = element_text(size = 10),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 45),
          strip.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.position = "bottom")


ggsave(bt_plt,
       filename = here::here("figs/fig1_bt_plt.svg"),
       dpi = 300,
       height = 8,
       width = 6)

## tidal current maximum----
tc_pal <- pals::brewer.prgn(n)
tc_plt <-
  ggplot() +
  geom_sf(data = agg_tc, aes(fill = tidal_curr, color = tidal_curr)) +
    scale_fill_gradientn(colors = tc_pal) +
    scale_color_gradientn(colors = tc_pal) +
    guides(fill = guide_colorbar(
      even.steps = TRUE,
      barheight = unit(0.25, "in"),
      barwidth  = unit(3, "in"),
      direction     = "horizontal",
      title.position = "left",
      title.vjust   = 0.8
    ),
    color = "none") +
    labs(fill = expression("Max tidal current (cm"~s^-1*")")) +
  ggnewscale::new_scale_color() +
  geom_sf(data = manage_sf,
          aes(color = name),
          fill = "transparent",
          linewidth = 0.5,
          show.legend = F) +
  scale_color_manual(values = c("black","purple"))+
  ggnewscale::new_scale_color() +
  geom_sf(data = nbbtca, color = "purple",
          linewidth = 0.5,
          linetype= 2,
          show.legend = F,
          fill = "transparent") +

    # land---
    geom_sf(data = ak_land) +
    # survey locations---
    geom_point(data = surv_locs,
               aes(y = lat, x = lon),
               shape = "+",
               size = 2,
               alpha = 0.75) +
    # tagging data---
  geom_point(data = male_crab_movement,
             aes(y = deploy_lat,
                 x = deploy_lon),
             size = 0.5,
             color = "grey10",
             alpha = 0.75,
             show.legend = F) +
  geom_segment(data = male_crab_movement,
               aes(y = deploy_lat, x = deploy_lon,
                   yend = lat0, xend = lon0),
               arrow = arrow(length = unit(0.05, "inches")), show.legend = F,
               color = "grey10",
               alpha = 0.75,
               linewidth = 0.25) +
  labs(y = "Latitude",
       x = "Longitude") +
  # themes etc---
  theme_bw() +
  coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
           ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
           expand = FALSE) +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(angle = 0))

ggsave(tc_plt,
       filename = here::here("figs/fig1_tc_plt.svg"),
       dpi = 300,
       height = 6,
       width = 6)


# subset area from globe
glob_crs <- "+proj=ortho +lon_0=-163 +lat_0=56.72769"

# background for the globe
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons
world <- rnaturalearthhires::countries10 %>%
  st_make_valid() %>%
  st_intersection(ocean %>% st_transform(4326)) %>%
  st_transform(crs = glob_crs )

polygon_sf_bound <- polygon_sf %>%
  st_transform(., crs = st_crs(ocean))

bbox_bound <- st_bbox(polygon_sf_bound)

x_range <- bbox_bound$xmax - bbox_bound$xmin
y_range <- bbox_bound$ymax - bbox_bound$ymin


xlim_vals <- c(bbox_bound$xmin - 1200000,
               bbox_bound$xmax + 900000)
ylim_vals <- c(bbox_bound$ymin - 300000,
               bbox_bound$ymax + 1300000)

outer_world <- ggplot() +
  geom_sf(data = ocean, fill = "#c9d7ebff", color = NA) +
  geom_sf(data = manage_sf %>%
            st_transform(., st_crs(ocean)),
          aes(color = name),
          fill = "transparent",
          linewidth = 0.2,
          show.legend = F) +
  scale_color_manual(values = c("black","purple")) +
    geom_sf(data = nbbtca, color = "purple",
            linewidth = 0.5,
            show.legend = F,linetype = "11",
            fill = "transparent") +
  geom_sf(data = world,linewidth = .2) +
  coord_sf(xlim = xlim_vals, ylim = ylim_vals) +
    annotation_scale(location = "bl") +
    annotation_north_arrow(which_north = "true",
                           style = north_arrow_minimal(),
                           location = "bl",
                           pad_y = unit(10, "mm")) +
  theme_fade() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(outer_world,
       file = here::here("figs/fig1_world.svg"),
       width = 6,
       height = 6)

