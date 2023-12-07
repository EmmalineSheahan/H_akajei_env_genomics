# script to produce project figures

library(ggplot2)
library(sf)
library(sp)
library(terra)
library(dplyr)
library(ggspatial)
library(paletteer)
library(RColorBrewer)
library(voluModel)
library(raster)

# create land for plotting, create target CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)

# reading in data
akajei_occs_down <- read.csv('./data/akajei_occs_down.csv')
site_metadata <- read.csv('./data/site_metadata.csv')
h_akajei_threshold <- rast('./results/h_akajei_threshold.tif')
h_akajei_suitability <- rast('./results/h_akajei_suitability.tif')
envs_for_plotting <- rast('./Env_data/envs_in_model_at_surface.tif')

# organizing site metadata
site_metadata1 <- site_metadata %>% select(c("lon", "lat", "H_akajei_presence"))
data_type <- rep("eDNA", times = nrow(site_metadata1))
site_metadata1 <- cbind(site_metadata1, data_type)
site_metadata1$H_akajei_presence[which(site_metadata1$H_akajei_presence 
                                      == 0)] <- "Absent"
site_metadata1$H_akajei_presence[which(site_metadata1$H_akajei_presence 
                                      == 1)] <- "Present"
colnames(site_metadata1) <- c("longitude", "latitude", "presence", "data_type")

# adding data type and presence to occurrences
data_type <- rep("Occurrence", times = nrow(akajei_occs_down))
presence <- rep("Present", times = nrow(akajei_occs_down))
akajei_occs_down_all <- cbind(akajei_occs_down, presence, data_type)

# merging data for plotting
full_df <- rbind(akajei_occs_down_all[,-3], site_metadata1)

# finding extent for larger map
akajei_sp <- akajei_occs_down
coordinates(akajei_sp) <- ~longitude+latitude
proj4string(akajei_sp) <- target_crs

# finding extent for inset map
akajei_small <- data.frame(site_metadata$lon, site_metadata$lat)
colnames(akajei_small) <- c("longitude", "latitude")
coordinates(akajei_small) <- ~longitude+latitude
proj4string(akajei_small) <- target_crs

# figure of occurrences
pdf('./plots/occurrence_and_eDNA_plot.pdf', width = 10, height = 5)
p <- ggplot(data = land) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5, pad_x = unit(2.55, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(4.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(119, 141), ylim = c(22, 36)) +
  geom_point(data = full_df, aes(x = longitude, y = latitude, col = presence, 
                                 shape = data_type), 
             size = 2.2) +
  scale_shape_manual(values=c(2, 16)) +
  scale_color_manual(values=c("#A84268", 'darkgreen')) +
  theme(panel.background = element_rect(fill = "lightblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(colour = "Presence", shape = "Data Origin")
print(p)
dev.off()

pdf('./plots/occurrence_and_edna_inset.pdf')
p <- ggplot(data = land) +
  geom_sf() +
  coord_sf(xlim = c(138, 141), ylim = c(34, 36)) +
  geom_point(data = full_df, aes(x = longitude, y = latitude, col = presence, 
                                 shape = data_type), 
             size = 2.2) +
  scale_shape_manual(values=c(2, 16)) +
  scale_color_manual(values=c("#A84268", 'darkgreen')) +
  theme(panel.background = element_rect(fill = "lightblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(colour = "Presence", shape = "Data Origin")
print(p)
dev.off()

# plots of environmental variables used for plotting
pdf('./plots/conductivity.pdf')
plot(envs_for_plotting[[1]], main = "Winter Conductivity (S/m)")
plot(land_poly, col = "grey", add = T)
dev.off()

pdf('./plots/salinity_summer.pdf')
plot(envs_for_plotting[[2]], main = "Summer Salinity (unitless)")
plot(land_poly, col = "grey", add = T)
dev.off()

pdf('./plots/salinity_winter.pdf')
plot(envs_for_plotting[[3]], main = "Winter Salinity (unitless)")
plot(land_poly, col = "grey", add = T)
dev.off()

pdf('./plots/winter_temperature.pdf')
plot(envs_for_plotting[[1]], main = "Winter Temperature (°C)")
plot(land_poly, col = "grey", add = T)
dev.off()

# layerplot of model for H akajei
pdf('./plots/h_akajei_levelplot.pdf')
plotLayers(rast = h_akajei_threshold, land = land, landCol = "grey", 
           title = "Predicted Suitable Habitat for H. akajei")
dev.off()

# plot of suitability and eDNA
only_present <- site_metadata1 %>% filter(presence == "Present")
only_absent <- site_metadata1 %>% filter(presence == "Absent")
coordinates(only_present) <- ~longitude+latitude
coordinates(only_absent) <- ~longitude+latitude
proj4string(only_present) <- target_crs
proj4string(only_absent) <- target_crs
ras_suit <- raster(h_akajei_suitability[[1]])
ras_suit <- crop(ras_suit, akajei_small)
ras_thresh <- raster(h_akajei_threshold[[1]])
ras_thresh <- crop(ras_thresh, akajei_small)

pdf('./plots/points_and_suitability.pdf')
plot(ras_suit)
plot(land_poly, col = "grey", add = T)
plot(only_absent, col = "#A84268", pch = 16, add = T)
plot(only_present, col = 'darkgreen', pch = 16, add = T)
dev.off()

pdf('./plots/points_and_threshold.pdf')
plot(ras_thresh)
plot(land_poly, col = "grey", add = T)
plot(only_absent, col = "#A84268", pch = 16, add = T)
plot(only_present, col = 'darkgreen', pch = 16, add = T)
dev.off()
