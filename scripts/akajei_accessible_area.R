# accessible area for H. akajei

library(sf)
library(sp)
library(raster)
library(terra)
library(voluModel)
library(rgdal)
library(dplyr)
library(rgeos)
library(rnaturalearthhires)
library(ENMeval)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4
cea_crs <- "+proj=cea +lat_ts=0 +lon_0"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)
land_poly_cea <- spTransform(land_poly, CRSobj = cea_crs)

coasts <- st_as_sf(rnaturalearthhires::coastline10)[1] %>% 
  st_set_crs(target_crs)
coast_poly <- as_Spatial(coasts)
coast_poly_cea <- spTransform(coast_poly, CRSobj = cea_crs)

# reading in downsampled occurrences
h_akajei <- read.csv('./data/akajei_occs_down.csv')

# reading in environmental data
temperature <- rast('./Env_data/temperature_summer_18.tif')
salinity <- rast('./Env_data/salinity_summer_18.tif')
conductivity <- rast('./Env_data/conductivity_summer_18.tif')
density <- rast('./Env_data/density_summer_18.tif')

# there are 7 occurrences for depths around 40, 45, and 50 m; enough to run a maxent
# model via ENMeval

wanted_env <- c(temperature[[10]], salinity[[10]], conductivity[[10]], 
                density[[10]])
for (i in 1:dim(wanted_env)[3]) {
  wanted_env[[i]] <- interpolateRaster(wanted_env[[i]])
}

# creating a depth list
wanted_depths <- unique(h_akajei$depth)

# creating various accessible areas of buffer distances 4, 6, 8, and 10 times 
# larger than the area of the coordinates, and then running initial models to 
# select the best buffer for each depth slice via AICc
buff_props <- c(4, 6, 8, 10)

# retrieving wanted occurrences
needed_slices <- c(wanted_depths[8], wanted_depths[9], wanted_depths[10])
new_df1 <- h_akajei %>% filter(depth == needed_slices[1])
new_df2 <- h_akajei %>% filter(depth == needed_slices[2])
new_df3 <- h_akajei %>% filter(depth == needed_slices[3])
new_df <- rbind(new_df1, new_df2, new_df3)
new_df_sp <- unique(new_df[,1:2])

# creating polygon and transforming to cea to deal in meters
temp_alph <- marineBackground(new_df_sp, fraction = 1, 
                              partCount = 2, 
                              clipToOcean = T,
                              buff = 600000)
temp_alph <- spTransform(as_Spatial(st_as_sf(temp_alph)), 
                         CRSobj = cea_crs)

# creating a list of shapes for each buffer proportion
shape_list <- vector("list", length = length(buff_props))
for (i in seq_along(buff_props)) {
  buffDist = (sqrt(buff_props[i]*gArea(temp_alph)) - sqrt(gArea(temp_alph)))/2
  shape_new <- raster::buffer(x = temp_alph, width = buffDist, dissolve = T)
  shape_new_clipped <- gDifference(shape_new, land_poly_cea)
  shape_sf <- st_as_sf(shape_new_clipped)
  shape_final <- st_transform(shape_sf, target_crs)
  shape_list[[i]] <- shape_final
}

# i loop to crop rasters to each polygon, run ENMeval for each polygon, and store
# model results for each polygon
model_results <- vector(length = length(shape_list))
for (i in 1:length(shape_list)) {
  cropped_env <- crop(wanted_env, shape_list[[i]])
  clipped_env <- mask(cropped_env, shape_list[[i]])
  simple_model <- ENMevaluate(occs = new_df_sp, envs = clipped_env, 
                              algorithm = "maxnet",
                              tune.args = list(fc = c("L"), rm = 1),
                              partitions = "none", n.bg = 10000)
  mod_res <- min(eval.results(simple_model)$AICc)
  model_results[i] <- mod_res
}

# find which polygon produced the lowest AIC
want_shape_num <- which(model_results == min(model_results, na.rm = T))

# for loop to create accessible areas per depth slice
depth_acc_area_list <- vector("list", length = length(wanted_depths))
for (i in seq_along(wanted_depths)) {
  
  # pull all occurrences for a given depth slice. for each depth slice, we want
  # occurrences from the slice above and the slice below in order to create 
  # connectivity
  if(i == 1) {
    needed_slices <- c(wanted_depths[i], wanted_depths[i+1])
    new_df1 <- h_akajei %>% filter(depth == needed_slices[1])
    new_df2 <- h_akajei %>% filter(depth == needed_slices[2])
    new_df <- rbind(new_df1, new_df2)
  } else if(i == length(wanted_depths)) {
    needed_slices <- c(wanted_depths[i-1], wanted_depths[i])
    new_df1 <- h_akajei %>% filter(depth == needed_slices[1])
    new_df2 <- h_akajei %>% filter(depth == needed_slices[2])
    new_df <- rbind(new_df1, new_df2)
  } else {
    needed_slices <- c(wanted_depths[i-1], wanted_depths[i], wanted_depths[i+1])
    new_df1 <- h_akajei %>% filter(depth == needed_slices[1])
    new_df2 <- h_akajei %>% filter(depth == needed_slices[2])
    new_df3 <- h_akajei %>% filter(depth == needed_slices[3])
    new_df <- rbind(new_df1, new_df2, new_df3)
  }
  new_df_sp <- unique(new_df[,1:2])
  
  # creating polygon and transforming to cea to deal in meters
  temp_alph <- marineBackground(new_df_sp, fraction = 1, 
                                partCount = 2, 
                                clipToOcean = T,
                                buff = 600000)
  temp_alph <- spTransform(as_Spatial(st_as_sf(temp_alph)), 
                           CRSobj = cea_crs)
  
  # buffering shape by a fact of 4 times the area covered by the points
  buffDist = (sqrt(buff_props[want_shape_num]*gArea(temp_alph)) - 
                sqrt(gArea(temp_alph)))/2
  shape_new <- raster::buffer(x = temp_alph, width = buffDist, dissolve = T)
  shape_new_clipped <- gDifference(shape_new, land_poly_cea)
  shape_sf <- st_as_sf(shape_new_clipped)
  shape_final <- st_transform(shape_sf, target_crs)
  depth_acc_area_list[[i]] <- shape_final
}

depth_acc_area <- do.call(rbind, depth_acc_area_list)

st_write(depth_acc_area, dsn = './Accessible_area', 
         layer = 'h_akajei_acc_area',
         driver = "ESRI Shapefile",
         append = F)


