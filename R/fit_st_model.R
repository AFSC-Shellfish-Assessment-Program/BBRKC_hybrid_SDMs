library(tidyverse)
library(raster)
library(patchwork)
library(sf)
library(stars)
library(sdmTMB)
library(INLA)
library(concaveman)

# Load data----
ak_crs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# boundary and grid data
load(here::here("data/spatial_layers.rdata"))

# load temperature data (october bottom temps from ROMS)----
load(here::here('data/agg_temp_interannual_sum_aut_mom6.rdata'))

# load bathymetry data----
load(here::here("data/agg_depth.rdata"))

# load sediment grain size data----
load(here::here("data/agg_phi.rdata"))

# load maxmimum tidal current data----
load(here::here("data/agg_tc.rdata"))

# load movement model grid for visualization
load(here::here("data/movement_model_particulars.rdata"))

# helper functions
source(here::here("R/helpers.R"))

# prediction grid (taken from merged environmental data to ensure it maps identically to
# movement model outputs)
pred_grid <-
  agg_temp_interannual_sum_aut %>%
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

# survey data-----
df_sf <-
  read_csv(file = here::here("data/RKC_matmale_matfem_2005-2024.csv")) %>%
  dplyr::rename_all(str_to_lower) %>%
  dplyr::rename(year = survey_year,
                lon = longitude,
                lat = latitude) %>%
  mutate(cpue = cpue/3.429) %>% # convert CPUE from nm^-2 to km^-2
  filter(year %in% 2005:2023,
         mat_sex == "Mature Male") %>%
  st_as_sf(., coords = c("lon", "lat"),
           crs = 4326) %>%
  st_transform(., crs = ak_crs)

# as data frame
df <- df_sf %>%
  sfc_as_cols(., names = c("lon","lat")) %>%
  st_set_geometry(NULL) %>%
  mutate(yf = factor(year))

# mesh creation----
survey_buffed <-
  df_sf %>%
  concaveman::concaveman(.) %>%
  st_union() %>%
  st_buffer(dist = 15) %>%
  as_Spatial()

max.edge <- 25
bound.outer <- 30

## INLA mesh----
mesh <- INLA::inla.mesh.2d(boundary = survey_buffed,
                           max.edge = c(1, 5) * max.edge,
                           offset = c(max.edge, bound.outer),
                           cutoff = max.edge)

## convert for sdmTMB----
ebs_spde <- sdmTMB::make_mesh(data = df,
                              xy_cols = c("lon", "lat"),
                              mesh = mesh)

## Visualizing for appendix------------
a <-   st_bbox(c(xmin = -1500, xmax = -50,
                 ymin = 400, ymax = 1500),
               crs = st_crs(ak_land)) %>%
  st_as_sfc() %>%
  st_make_grid(., cellsize = 150)

# mesh figure-----
png(filename = here::here("figs/ebs_mesh.png"),
    width = 10, height = 8,
    units = "in", res = 300)

plot(a,
     col = "white",
     axes = T,
     xlim = c(-1000, -500),
     ylim = c(500, 1550),
     xlab = "Eastings (km)",
     ylab =  "Northings (km)",
     xaxt = "n",
     yaxt = "n",
     border = "grey90",
     cex.lab = 1.7,
     cex.axis = 1.5)
# plot(ebs_grid, add = T)
plot(ak_land,
     col = "grey",
     xlim = c(-1000, -500),
     ylim = c(500, 1400),
     add = T)
plot(ebs_spde, add = T)
axis(1, at = seq(-1500, 50, by = 150))
axis(3,at = seq(-1500, 50, by = 150), labels = F, tick = T)
axis(2, at = seq(400, 1500, by = 150))
axis(4, labels = F, tick = T,at = seq(400, 1500, by = 150))
dev.off()

# prediction grid for spatiotemporal and movement models----
pred_grid_sfc <- pred_grid %>%
  filter(year == 2023, month == 10) %>%
  st_as_sfc()

# prediction grid figure----
png(filename = here::here("figs/pred_grid.png"),
    width = 10, height = 8,
    units = "in", res = 300)
plot(a,
     col = "white",
     axes = T,
     xlim = c(-1000, -500),
     ylim = c(500, 1400),
     xlab = "Eastings (km)",
     ylab =  "Northings (km)",
     xaxt = "n",
     yaxt = "n",
     border = "grey90",
     cex.lab = 1.7,
     cex.axis = 1.5)
plot(ak_land,
     col = "grey",
     xlim = c(-1000, -500),
     ylim = c(500, 1400),
     add = T)
plot(pred_grid_sfc, add = T, border = "purple", lwd = 3)
plot(grid2, add = T, border = "black")
axis(1, at = seq(-1500, 50, by = 150))
axis(3,at = seq(-1500, 50, by = 150), labels = F, tick = T)
axis(2, at = seq(400, 1500, by = 150))
axis(4, labels = F, tick = T,at = seq(400, 1500, by = 150))
dev.off()

## fit IID tweedie glmm w/ spatial and spatiotemporal random fields----
m1_sdm <- sdmTMB::sdmTMB(
  formula = cpue ~ 0 + as.factor(year),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  mesh = ebs_spde,
  family = sdmTMB::tweedie(),
  data = df
)

summary(m1_sdm)
sanity(m1_sdm)

# prediction grid by year-----
nd <- pred_grid %>%
  st_centroid() %>%
  dplyr::select(year) %>%
  sfc_as_cols(., names = c("lon","lat")) %>%
  st_set_geometry(NULL) %>%
  filter(year != 2020) # no survey in 2020

# make predictions--------
pred <- sdmTMB:::predict.sdmTMB(object = m1_sdm,
                                newdata = nd,
                                type = "response") %>%
  st_as_sf(., coords = c("lon","lat"), crs = ak_crs)

# predictions as sf----
pred_df_sf <-
  pred_grid %>%
  dplyr::select(geometry) %>%
  st_join(., pred) #%>%

df$resids <- residuals(m1_sdm,
                       type = "mle-mvn")
ggplot(df) +
  geom_histogram(aes(resids)) +
  facet_wrap(~year) +
  labs(title = "Spatiotemporal model residuals by year") +
  theme_minimal()

# save outputs----
save(pred,
     df,
     pred_df_sf,
     m1_sdm,
     pred_grid,
     nd,
     file = here::here("data/st_model_predictions.rdata"))
