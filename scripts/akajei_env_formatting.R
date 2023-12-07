# formatting environmental variables for H akajei project

# preparing environmental data for 3D and 2D modeling

library(sf)
library(terra)
library(dplyr)
library(voluModel)
library(R.utils)
library(rnaturalearth)
library(rnaturalearthhires)
library(rgdal)
library(ncdf4)
library(raster)
library(rgeos)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4
cea_crs <- "+proj=cea +lat_ts=0 +lon_0"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)
land_poly_cea <- spTransform(land_poly, CRSobj = cea_crs)

# reading in occurrences for extent matching
h_akajei <- read.csv('./data/akajei_occs_down.csv')

# creating base raster for standardized extent, cropped to the extent 
# of a shape buffered for all of the points
temp_alph <- marineBackground(h_akajei[,1:2], fraction = 1, 
                              partCount = 2, 
                              clipToOcean = T,
                              buff = 600000)
temp_alph <- spTransform(as_Spatial(st_as_sf(temp_alph)), 
                         CRSobj = cea_crs)

# buffering shape by a fact of 6 times the area covered by the points
buffDist <- (sqrt(8*gArea(temp_alph)) - sqrt(gArea(temp_alph)))/2
shape_new <- raster::buffer(x = temp_alph, width = buffDist, dissolve = T)
shape_new_clipped <- gDifference(shape_new, land_poly_cea)
shape_sf <- st_as_sf(shape_new_clipped)
shape_final <- st_transform(shape_sf, target_crs)
wanted_extent <- ext(shape_final)

temp_csv <- read.csv('./Env_data/woa18_decav_t13mn04.csv', header = F)
depth_vals <- c(0, matrix(unlist(temp_csv[2, 4:ncol(temp_csv)]))[,1])
for (i in 1:length(depth_vals)) {
  depth_vals[i] <- paste0("X", depth_vals[i])
}
temp_csv <- temp_csv[-(1:2),]
new_colnames <- c("Latitude", "Longitude", depth_vals)
colnames(temp_csv) <- new_colnames
temp_df <- temp_csv %>% dplyr::select(Longitude, Latitude, depth_vals[11])
temp_rast_ext <- rast(temp_df, type = "xyz", crs = target_crs, 
                      extent = wanted_extent)

# env_maker function
# wanted_data = filepath of csv with env data
# output_name = filepath of output tif
# converts WOA csv to spatraster stack of depthslices
env_maker <- function(wanted_data, output_name) {
  temp_csv <- read.csv(wanted_data, header = F)
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
    temp_rast <- rast(temp_df, type = "xyz", crs = target_crs, 
                      extent = ext(temp_rast_ext))
    rast_list[[i]] <- temp_rast
  }
  final_rast <- rast(rast_list)
  # I'm buffering the lowest occurrence depth by 100 m
  final_rast <- final_rast[[c(1:27)]]
  writeRaster(final_rast, filename = output_name, overwrite = T)
}

# Temperature
env_maker(wanted_data = './Env_data/woa18_decav_t13mn04.csv', 
          output_name = './Env_data/temperature_winter_18.tif')
env_maker(wanted_data = './Env_data/woa18_decav_t15mn04.csv',
          output_name = './Env_data/temperature_summer_18.tif')

# Salinity
env_maker(wanted_data = './Env_data/woa18_decav_s13mn04.csv', 
          output_name = './Env_data/salinity_winter_18.tif')
env_maker(wanted_data = './Env_data/woa18_decav_s15mn04.csv',
          output_name = './Env_data/salinity_summer_18.tif')

# Density
env_maker(wanted_data = './Env_data/woa18_decav_I13mn04.csv', 
          output_name = './Env_data/density_winter_18.tif')
env_maker(wanted_data = './Env_data/woa18_decav_I15mn04.csv',
          output_name = './Env_data/density_summer_18.tif')

# Conductivity
env_maker(wanted_data = './Env_data/woa18_A5B7_C13mn04.csv', 
          output_name = './Env_data/conductivity_winter_18.tif')
env_maker(wanted_data = './Env_data/woa18_A5B7_C15mn04.csv',
          output_name = './Env_data/conductivity_summer_18.tif')