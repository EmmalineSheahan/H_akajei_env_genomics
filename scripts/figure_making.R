# script to produce project figures

library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(ggspatial)
library(paletteer)

# create land for plotting, create target CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)

# reading in data
akajei_occs_down <- read.csv('./data/akajei_occs_down.csv')
site_metadata <- read.csv('./data/site_metadata.csv')

# organizing site metadata
site_metadata <- site_metadata %>% select(c("lon", "lat", "H_akajei_presence"))
data_type <- rep("eDNA", times = nrow(site_metadata))
site_metadata <- cbind(site_metadata, data_type)
site_metadata$H_akajei_presence[which(site_metadata$H_akajei_presence 
                                      == 0)] <- "Absent"
site_metadata$H_akajei_presence[which(site_metadata$H_akajei_presence 
                                      == 1)] <- "Present"
colnames(site_metadata) <- c("longitude", "latitude", "presence", "data_type")

# adding data type and presence to occurrences
data_type <- rep("Occurrence", times = nrow(akajei_occs_down))
presence <- rep("Present", times = nrow(akajei_occs_down))
akajei_occs_down <- cbind(akajei_occs_down, presence, data_type)


# merging data for plotting
full_df <- rbind(akajei_occs_down[,-3], site_metadata)

# finding extent
akajei_sp <- akajei_occs_down
coordinates(akajei_sp) <- ~longitude+latitude
proj4string(akajei_sp) <- target_crs

# figure of occurrences
pdf('./plots/occurrence_and_eDNA_plot.pdf')
ggplot(data = land) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5, pad_x = unit(3.55, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(6.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(119, 140), ylim = c(22, 36)) +
  geom_point(data = full_df, aes(x = longitude, y = latitude, col = presence, 
                                 shape = data_type), 
             size = 2.2) +
  scale_shape_manual(values=c(2, 16)) +
  scale_color_manual(values=c("#A84268", 'darkgreen')) +
  theme(panel.background = element_rect(fill = "lightblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(colour = "Presence", shape = "Data Origin")
dev.off()
