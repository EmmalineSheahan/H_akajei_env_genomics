# 3D ecological niche modelling for H akajei

source('./scripts/all_functions_3D_sdm.R')

library(raster)
library(dplyr)
library(maxnet)
library(rnaturalearthhires)
library(dismo)
library(ENMeval)
library(voluModel)
library(sf)
library(usdm)
library(doParallel)

# create land for plotting, create target CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)

# reading in occurrences
h_akajei <- read.csv('./data/akajei_occs_down.csv')

# reading in environmental rasters
density_summer_18 <- rast('./Env_data/density_summer_18.tif')
density_winter_18 <- rast('./Env_data/density_winter_18.tif')
conductivity_summer_18 <- rast('./Env_data/conductivity_summer_18.tif')
conductivity_winter_18 <- rast('./Env_data/conductivity_winter_18.tif')
salinity_summer_18 <- rast('./Env_data/salinity_summer_18.tif')
salinity_winter_18 <- rast('./Env_data/salinity_winter_18.tif')
temperature_summer_18 <- rast('./Env_data/temperature_summer_18.tif')
temperature_winter_18 <- rast('./Env_data/temperature_winter_18.tif')

envs_list <- list(density_summer_18, density_winter_18, conductivity_summer_18,
                  conductivity_winter_18, salinity_summer_18,
                  salinity_winter_18, temperature_summer_18, 
                  temperature_winter_18)
envs_names <- c("density_summer", "density_winter", "conductivity_summer", 
                "conductivity_winter", "salinity_summer",
                "salinity_winter", "temperature_summer", "temperature_winter")

# read in accessible area 
acc_area <- st_read(dsn = './Accessible_area', layer = 'h_akajei_acc_area')

# converting envs_list into list of stacks where each element is a depth layer
envs_new <- env_stack_transform(envs_list, envs_names)

# cropping environmental variables to the accessible area on a per depth slice
# basis

# creating depth slice list to match between envs and acc areas
depths <- unique(h_akajei$depth)
env_depths <- as.numeric(gsub("X", '', names(envs_list[[1]])))
acc_area_call <- vector("list", length = length(depths))
for(i in 1:length(depths)) {
  if(i == 1) {
    t <- env_depths[which(env_depths <= depths[i])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  } else if(i == length(depths)) {
    t <- env_depths[which(env_depths > depths[i-1])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  } else {
    t <- env_depths[which(env_depths <= depths[i] & env_depths > depths[i-1])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  }
}
acc_area_call <- unlist(acc_area_call)

# using crop_stack function to mask each depth layer to the associated accessible
# area
envs_final <- vector("list", length = length(envs_new))
for(i in 1:length(envs_new)) {
  envs_final[[i]] <- crop_stack(enviro = envs_new[[i]], 
                                accarea = acc_area[acc_area_call[i],], 
                                which_interpolate = c(1:8))
}

# creating list of spatial points by depth slice
h_akajei_sp_list <- xyzmat_to_3Dsp(depth_slices = env_depths,
                                      occ_mat = h_akajei)

# creating list of background points by depth slice
h_akajei_bg_list <- bg_list_maker(enviro_stack = envs_final,
                                     depth_slices = env_depths,
                                     wanted_var = "temperature_summer",
                                     wanted_num = 1000,
                                     presences = h_akajei_sp_list)

# selecting variables with highest permutation importance and which reduce VIF
envs_for_model <- var_select(wanted_brick = envs_final,
                             wanted_sp = h_akajei_sp_list,
                             wanted_bg = h_akajei_bg_list)
writeRaster(envs_for_model[[1]], filename = './Env_data/envs_in_model_at_surface.tif')

# producing 3D models
# extracting values at occs and bg
all_h_akajei_df <- brick_extract(wanted_brick = envs_for_model, 
                                    wanted_sp = h_akajei_sp_list,
                                    wanted_bg = h_akajei_bg_list)

# defining partition schemes
h_akajei_k <- partition_3D(h_akajei_sp_list, h_akajei_bg_list, 
                              all_h_akajei_df[[3]],
                              'k.fold',
                              kfolds = 5)

h_akajei_block <- partition_3D(h_akajei_sp_list, h_akajei_bg_list, 
                                  all_h_akajei_df[[3]],
                                  'block',
                                  orientation = 'lat_lon')

# training and testing models
h_akajei_mod_k <- maxent_3D(df_for_maxent = all_h_akajei_df[[3]],
                               wanted_fc = c("L", "LQ", "LQH", "LQP"),
                               wanted_rm = c(1:4),
                               wanted_partition = h_akajei_k,
                               projection_layers = envs_for_model,
                               occ_list = h_akajei_sp_list)
writeRaster(h_akajei_mod_k$predictions[[12]], 
            filename = './results/h_akajei_suitability.tif')
write.csv(h_akajei_mod_k$results, file = './results/h_akajei_mod_results.csv')

# LQH with rm 4 has the lowest AICc

# thresholding
h_akajei_threshold <- threshold_3D(projection_layers = NA, 
                                   predicted_layers = h_akajei_mod_k$predictions[[12]],
                                   thresholding_vals = c(0.90, 0.95), 
                                   occ_list = h_akajei_sp_list,
                                   bg_list = h_akajei_bg_list)
writeRaster(h_akajei_threshold$threshold_layers, 
            filename = './results/h_akajei_threshold.tif', overwrite = T)
write.csv(h_akajei_threshold$tss_results, file = './results/h_akajei_tss.csv')

# running a model evaluation with eDNA test points and absences
site_metadata <- read.csv('./data/site_metadata.csv')
wanted_presence <- site_metadata %>% filter(H_akajei_presence == 1)
wanted_absence <- site_metadata %>% filter(H_akajei_presence == 0)
wanted_presence <- data.frame(wanted_presence$lon, wanted_presence$lat)
wanted_absence <- data.frame(wanted_absence$lon, wanted_absence$lat)
eDNA_pres <- terra::extract(envs_for_model[[1]], wanted_presence)[,-1]
eDNA_abs <- terra::extract(envs_for_model[[1]], wanted_absence)[,-1]
#psuedo_abs <- all_h_akajei_df[[3]] %>% filter(presence == 0)
#eDNA_abs <- rbind(eDNA_abs, psuedo_abs[,-1])
eDNA_eval <- dismo::evaluate(p = eDNA_pres, a = eDNA_abs, 
                             model = h_akajei_mod_k$models[[12]])

# calculate sensitivity, specificity, and tss for eDNA sites
confirmed_pres <- terra::extract(h_akajei_threshold$threshold_layers[[1]], 
                                 wanted_presence)
sensitivity <- (length(which(confirmed_pres$mean == 1))/nrow(wanted_presence))
confirmed_abs <- terra::extract(h_akajei_threshold$threshold_layers[[1]], 
                                wanted_absence)
specificty <- (length(which(confirmed_abs$mean == 0))/nrow(wanted_absence))
eDNA_tss <- (sensitivity + ((1/3)*specificity)) - 1
