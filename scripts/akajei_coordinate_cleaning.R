# gathering and cleaning occurrence records for Hemitrygon akajei

library(dismo)
library(robis)
library(rnaturalearthhires)
library(dplyr)
library(sf)
library(terra)
library(voluModel)
library(CoordinateCleaner)

# create land for plotting, create target CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)

# pulling in raw data from GBIF
akajei_gbif_raw <- gbif(genus = 'Hemitrygon', species = 'akajei', geo = T,
                           removeZeros = T, download = T)

# selecting only records with depth information
akajei_gbif_raw <- akajei_gbif_raw[complete.cases(akajei_gbif_raw$depth),]
akajei_gbif_raw <- akajei_gbif_raw[akajei_gbif_raw$depth >= 0.0,]
akajei_gbif_raw <- akajei_gbif_raw[akajei_gbif_raw$depth < 4000.0,]

# removing any duplicates
akajei_gbif_raw <- akajei_gbif_raw %>% select(lat, lon, depth, eventDate)
akajei_gbif <- distinct(akajei_gbif_raw)

# pulling in raw data from OBIS
akajei_obis_raw <- occurrence(scientificname = "Hemitrygon akajei", qcfields = T)

# only keeping data with locality information
akajei_obis_raw <- akajei_obis_raw[complete.cases(akajei_obis_raw$decimalLatitude),]
akajei_obis_raw <- 
  akajei_obis_raw[complete.cases(akajei_obis_raw$decimalLongitude),]

# only keeping data with depth information
akajei_obis_raw <- akajei_obis_raw[complete.cases(akajei_obis_raw$depth),]
akajei_obis_raw <- akajei_obis_raw[akajei_obis_raw$depth >= 0.0,]
akajei_obis_raw <- akajei_obis_raw[akajei_obis_raw$depth < 4000.0,]

# removing duplicates
akajei_obis_raw <- akajei_obis_raw %>% select(decimalLatitude, decimalLongitude, 
                                              depth, eventDate)
akajei_obis <- distinct(akajei_obis_raw)

# combining datasets and performing final cleaning
colnames(akajei_obis) <- colnames(akajei_gbif)
akajei_obis <- data.frame(akajei_obis)
akajei_occs <- rbind(akajei_gbif, akajei_obis)
akajei_occs <- distinct(akajei_occs)
akajei_occs <- akajei_occs[-which(is.na(akajei_occs$lat)),]
akajei_occs <- akajei_occs[,1:3]
colnames(akajei_occs) <- c("latitude", "longitude", "depth")

# removing outliers
species <- rep("Hemitrygon_akajei", times = nrow(akajei_occs))
akajei_occs <- cbind(akajei_occs, species)
akajei_occs_outl <- cc_outl(akajei_occs, lon = "longitude", 
                        lat = "latitude", species = "species", 
                        method = "distance", value = "flagged", 
                        tdi = 500)
akajei_occs_new <- akajei_occs[which(akajei_occs_outl),]
akajei_occs_sp <- akajei_occs_new
coordinates(akajei_occs_sp) <- ~longitude+latitude
proj4string(akajei_occs_sp) <- target_crs
akajei_occs_sf <- st_as_sf(akajei_occs_sp)

plot(akajei_occs_sp, col = "red")
plot(land_poly, col = NA, add = T)

# thinning records
# creating a raster with desired 1/4 degree resolution for thinning
temp_csv <- read.csv('./Env_data/woa18_decav_t15mn04.csv', header = F)
depth_vals <- c(0, matrix(unlist(temp_csv[2, 4:ncol(temp_csv)]))[,1])
for (i in 1:length(depth_vals)) {
  depth_vals[i] <- paste0("X", depth_vals[i])
}
temp_csv <- temp_csv[-(1:2),]
new_colnames <- c("Latitude", "Longitude", depth_vals)
colnames(temp_csv) <- new_colnames
rast_list <- vector("list", length = length(depth_vals))
for (i in 1:length(depth_vals)) {
  temp_df <- temp_csv %>% dplyr::select(Longitude, Latitude, depth_vals[i])
  temp_rast <- rast(temp_df, type = "xyz", crs = target_crs)
  rast_list[[i]] <- temp_rast
}
final_rast <- rast(rast_list)

# thinning occurrences
# Get the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("X", "", names(final_rast)))
akajei_occs_new$index <- unlist(lapply(akajei_occs_new$depth, 
                              FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(akajei_occs_new$index)

# downsampling occurrences
akajei_occs_down <- data.frame()
for(i in indices){
  tempPoints <- akajei_occs_new[akajei_occs_new$index==i,]
  tempPoints <- downsample(tempPoints, final_rast[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  akajei_occs_down <- rbind(akajei_occs_down, tempPoints)
}

# ordering by the depth column
akajei_occs_down <- akajei_occs_down[order(akajei_occs_down$depth, decreasing = F),]

# writing final occurrence records to file
write.csv(akajei_occs_down, file = './data/akajei_occs_down.csv', 
          row.names = F)
