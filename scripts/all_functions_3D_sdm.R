# all functions necessary for the 3D niche modeling pipeline
library(raster)
library(terra)
library(dplyr)
library(maxnet)
library(rnaturalearth)
library(dismo)
library(ENMeval)
library(voluModel)
library(sf)

# crop environmental layers to accessible area
# crop_stack function to crop all layers in stack
# enviro = environmental layer SpatRaster. each layer is a different 
# environmental variable
# which_interpolate = the numbers of layers to be interpolated
# acc_area = accessible area sf
# returns a raster stack cropped and masked at each depth layer to the 
# accessible area
crop_stack <- function(enviro, accarea, which_interpolate) {
    for (i in seq_along(which_interpolate)) {
      enviro[[i]] <- interpolateRaster(enviro[[i]])
    }
    temp <- mask(enviro, accarea)
  return(temp)
}

# transform occurrence into list of spatial points by depth slice
# xyzmat_to_3Dsp function to transform occurrence matrix into a list of 
# matrices with longitude and latitude headers by depth slice
# depth_slices = numeric vector of depth slices
# occ_mat = xyz occurrence matrix. expects x and y to be named "longitude" and 
# "latitude", and z to be named "depth."
xyzmat_to_3Dsp <- function(depth_slices, occ_mat) {
  sp_list <- vector("list", length = length(depth_slices))
  for (i in seq_along(depth_slices)) {
    test <- occ_mat %>% filter(depth >= depth_slices[i])
    test2 <- test %>% filter(depth < depth_slices[i] + 1)
    if (nrow(test2) == 0) {
      test3 <- NA 
    } else {
      test3 <- data.frame(test2$longitude, test2$latitude)
      colnames(test3) <- c("longitude", "latitude")
    }
    sp_list[[i]] <- test3
  }
  return(sp_list)
}

# transforming list of environmental raster stacks into a list of combined 
# stacks where each list element is a different depth layer. this is the
# input for bg_list_maker and brick_extract
# envs_all = list of SpatRaster stacks where each element is a different 
# environmental variable
env_stack_transform <- function(envs_all, envs_names) {
  env_call <- vector(length = length(envs_all))
  for (i in 1:length(envs_all)) {
    place <- paste0("envs_all[[", i, "]][[i]]")
    env_call[i] <- place
  }
  
  stack_list_all <- vector("list", length = dim(envs_all[[1]])[3])
  for (i in 1:dim(envs_all[[1]])[3]) {
    j_list <- vector("list", length = length(env_call))
    for (j in 1:length(env_call)) {
      j_list[[j]] <- eval(parse(text = env_call[j]))
    }
    st <- rast(j_list)
    names(st) <- envs_names
    stack_list_all[[i]] <- st
  }
  return(stack_list_all)
}

# creating a list of background points at each depth slice
# bg_list_maker function to create background per depth slice
# enviro_stack = list of SpatRaster stacks. each list element is a depth
# depth_slices = numeric vector of depth slice intervals
# wanted_var = name of the variable layer from which points will be based on
# wanted_num = desired number of background points
# presences = list of presence points by depth slice
# returns a list of spatial points at each depth slice
bg_list_maker <- function(enviro_stack, depth_slices, wanted_var, wanted_num, 
                          presences) {
  bg_list <- vector("list", length = length(depth_slices))
  for (i in 1:length(depth_slices)) {
    pull_this <- which(names(enviro_stack[[i]]) == wanted_var)
    bg <- dismo::randomPoints(raster(enviro_stack[[i]][[pull_this]]), 
                              n = wanted_num, p = presences[[i]]) %>% 
      as.data.frame()
    colnames(bg) <- c("longitude", "latitude")
    bg_list[[i]] <- bg
  }
  return(bg_list)
}

# Extract the environmental variables at the occurrence or 
# background points for each depth value
# brick_extract function outputs a list containing [[1]] an occs dataframe
# for ENMevaluate containing the longitude and latitudes of occurrences and 
# environmnetal values, [[2]] a bg dataframe for ENMevaluate containing the 
# longitudes and latitudes for background points and the associated environmental
# values, and [[3]] an environmental variable matrix with 
# presence values for maxent
# wanted_brick = list of SpatRaster stacks of environmental variables where each 
# element corresponds to a depth slice
# wanted_sp = list of occurrence coordinates by depth slice
# wanted_bg = list of background points by depth slice
brick_extract <- function(wanted_brick, wanted_sp, wanted_bg) {
  sp_part <- vector("list", length = length(wanted_brick))
  bg_part <- vector("list", length = length(wanted_brick))
  for (i in 1:length(wanted_brick)) {
    test_list_sp <- vector("list", length = length(names(wanted_brick[[i]])))
    test_list_bg <- vector("list", length = length(names(wanted_brick[[i]])))
    for(j in 1:length(names(wanted_brick[[i]]))) {
      if(any(is.na(wanted_sp[[i]]))) {
        test_sp <- matrix(nrow = 1, ncol = 1)
        test_sp[1,1] <- "delete.me"
        test_sp <- data.frame(test_sp)
        colnames(test_sp) <- names(wanted_brick[[i]][[j]])
        test_list_sp[[j]] <- test_sp
      } else {
      test_sp <- terra::extract(x = wanted_brick[[i]][[j]], y = wanted_sp[[i]])
      test_list_sp[[j]] <- test_sp %>% 
        dplyr::select(names(wanted_brick[[i]][[j]]))
      }
      test_bg <- terra::extract(x = wanted_brick[[i]][[j]], y = wanted_bg[[i]])
      test_list_bg[[j]] <- test_bg %>% 
        dplyr::select(names(wanted_brick[[i]][[j]]))
    }
    for(k in 1:length(wanted_sp)) {
      if(any(is.na(wanted_sp[[k]]))) {
        longitude <- NA
        latitude <- NA
        newdf <- data.frame(longitude, latitude)
        wanted_sp[[k]] <- newdf
      }
    }
    all_sp <- do.call(cbind, test_list_sp)
    all_sp <- cbind(wanted_sp[[i]], all_sp)
    all_bg <- do.call(cbind, test_list_bg)
    all_bg <- cbind(wanted_bg[[i]], all_bg)
    presence <- rep(1, times = nrow(all_sp))
    all_sp <- cbind(presence, all_sp)
    sp_part[[i]] <- all_sp
    presence <- rep(0, times = nrow(all_bg))
    all_bg <- cbind(presence, all_bg)
    bg_part[[i]] <- all_bg
  }
  final_sp <- do.call(rbind, sp_part)
  final_sp <- final_sp[-which(final_sp[,4] == "delete.me"),]
  enm_occs <- final_sp[,-1]
  final_bg <- do.call(rbind, bg_part)
  enm_bg <- final_bg[,-1]
  final_df <- rbind(final_sp[,-(2:3)], final_bg[,-(2:3)])
  all_products <- list(enm_occs, enm_bg, final_df)
  return(all_products)
}

# var_remove function for single iterative variable removal
# wanted_envs = list of spatraster stacks of environmental variables
# wanted_sp = list of occurrences by depth
# wanted_bg = list of background points by depth
# returns a list of spatrasters where the variable with the lowest permutation
# importance has been removed
var_remove <- function(wanted_envs, wanted_sp, wanted_bg) {
  maxent_df <- brick_extract(wanted_brick = wanted_envs, wanted_sp = wanted_sp,
                             wanted_bg = wanted_bg)
  initial_mod <- maxent(x = maxent_df[[3]][,-1], p = maxent_df[[3]][,1])
  initial_results <- data.frame(initial_mod@results)
  measure_name <- row.names(initial_results)
  perm_import <- data.frame(measure_name, 
                            initial_results[,1])
  colnames(perm_import) <- c("measure_name", "result")
  want_pull <- grep('.permutation.importance', perm_import$measure_name)
  final_perm_import <- perm_import[want_pull,]
  perm_remove <- which(final_perm_import$result == min(final_perm_import$result))
  if (length(perm_remove > 1)) {
    perm_remove <- perm_remove[1]
  }
  remove_this <- final_perm_import$measure_name[perm_remove]
  remove_this <- gsub('.permutation.importance', '', remove_this)
  
  # removing variable from analysis
  new_envs <- vector("list", length = length(wanted_envs))
  for(i in 1:length(wanted_envs)) {
    remove_env <- which(names(wanted_envs[[i]]) %in% remove_this) 
    new_envs[[i]] <- wanted_envs[[i]][[-remove_env]]
  }
  return(new_envs)
}

# var_select function to iterate var_remove until all vifs are below 5
# wanted_brick = list of environmental spatrasters
# wanted_sp = list of occurrences
# wanted_bg = list of background points
var_select <- function(wanted_brick, wanted_sp, wanted_bg) {
  all_depths <- vector("list", length = length(wanted_brick))
  for(i in 1:length(wanted_brick)) {
    all_envvars <- vector("list", length = length(names(wanted_brick[[i]])))
    for(j in 1:length(names(wanted_brick[[i]]))) {
      all_envvars[[j]] <- values(wanted_brick[[i]][[j]])[,1]
    }
    perdepth_df <- do.call(cbind, all_envvars)
    colnames(perdepth_df) <- names(wanted_brick[[i]])
    all_depths[[i]] <- perdepth_df
  }
  alldepth_env <- do.call(rbind, all_depths) %>% data.frame()
  t <- max(usdm::vif(alldepth_env)$VIF)
  final_env <- wanted_brick
  
  while(t > 10) {
    final_env <- var_remove(wanted_envs = final_env, wanted_sp = wanted_sp,
                            wanted_bg = wanted_bg)
    if(dim(final_env[[1]])[3] < 5) {
      t <- 4
    } else {
      all_depths <- vector("list", length = length(final_env))
      for(i in 1:length(final_env)) {
        all_envvars <- vector("list", length = length(names(final_env[[i]])))
        for(j in 1:length(names(final_env[[i]]))) {
          all_envvars[[j]] <- values(final_env[[i]][[j]])[,1]
        }
        perdepth_df <- do.call(cbind, all_envvars)
        colnames(perdepth_df) <- names(final_env[[i]])
        all_depths[[i]] <- perdepth_df
      }
      alldepth_env <- do.call(rbind, all_depths) %>% data.frame()
      t <- max(usdm::vif(alldepth_env)$VIF)
    }
  }
  
  return(final_env)
}

# partition_3D function creates 3D partitions given a specified partitioning scheme
# wanted_sp = list of coordinates of occurrence by depth
# wanted_bg = list of coordinates of background by depth
# maxent_df = [[3]] output of brick_extract
# which_partition = desired partitioning scheme as character
# orientation = for block. either 'lat_lon' or 'lon_lat'
# aggregation.factor = for checkerboard. what squares will be aggregated by
# kfolds = for kfold. desired number of folds
partition_3D <- function(wanted_sp, wanted_bg, maxent_df, which_partition, 
                         orientation,
                         aggregation.factor, kfolds) {
  partition_list <- vector("list", length = length(wanted_sp))
  if(which_partition == 'k.fold') {
    for(i in 1:length(wanted_sp)) {
      if(any(is.na(wanted_sp[[i]]))) {
        katdepth <- list(occs.grp = NA, bg.grp = sample(c(1:kfolds), 
                                                        nrow(wanted_bg[[i]]), 
                                                        replace = T))
        partition_list[[i]] <- katdepth
      } else {
        occs.grp <- sample(c(1:kfolds), nrow(wanted_sp[[i]]), replace = T)
        bg.grp <- sample(c(1:kfolds), nrow(wanted_bg[[i]]), replace = T)
        katdepth <- list(occs.grp = occs.grp, bg.grp = bg.grp)
        partition_list[[i]] <- katdepth
      }
    }
  } else if(which_partition == 'block') {
    for(i in 1:length(wanted_sp)) {
      if(any(is.na(wanted_sp[[i]]))) {
        blockatdepth <- get.block(occs = wanted_bg[[i]], bg = wanted_bg[[i]],
                                  orientation = orientation)
        blockatdepth$occs.grp <- NA
        partition_list[[i]] <- blockatdepth
      } else {
        blockatdepth <- get.block(occs = wanted_sp[[i]], bg = wanted_bg[[i]],
                                  orientation = orientation)
        partition_list[[i]] <- blockatdepth
      }
    }
  } else if(which_partition == 'checkerboard1') {
      for(i in 1:length(wanted_sp)) {
        if(any(is.na(wanted_sp[[i]]))) {
          checkatdepth <- get.checkerboard1(occs = wanted_bg[[i]], bg = wanted_bg[[i]],
                                    aggregation.factor = aggregation.factor)
          checkatdepth$occs.grp <- NA
          partition_list[[i]] <- checkatdepth
        } else {
          checkatdepth <- get.checkerboard1(occs = wanted_sp[[i]], 
                                          bg = wanted_bg[[i]],
                                  aggregation.factor = aggregation.factor)
          partition_list[[i]] <- checkatdepth
        }
      }
  } else if(which_partition == 'checkerboard2') {
    for(i in 1:length(wanted_sp)) {
      if(any(is.na(wanted_sp[[i]]))) {
        checkatdepth <- get.checkerboard2(occs = wanted_bg[[i]], bg = wanted_bg[[i]],
                                          aggregation.factor = aggregation.factor)
        checkatdepth$occs.grp <- NA
        partition_list[[i]] <- checkatdepth
      } else {
        checkatdepth <- get.checkerboard2(occs = wanted_sp[[i]], 
                                        bg = wanted_bg[[i]],
                                aggregation.factor = aggregation.factor)
        partition_list[[i]] <- checkatdepth
      }
    }
  } else {
    print("Not a usable partition")
  }
  occs_list <- vector("list", length = length(partition_list))
  bg_list <- vector("list", length = length(partition_list))
  for(i in 1:length(partition_list)) {
    occs_list[[i]] <- partition_list[[i]]$occs.grp
    bg_list[[i]] <- partition_list[[i]]$bg.grp
  }
  occ_partitions <- unlist(occs_list)[-which(is.na(unlist(occs_list)))]
  needed_presence <- maxent_df %>% filter(presence == 1)
  occ_partitions <- occ_partitions[which(complete.cases(needed_presence))]
  bg_partitions <- unlist(bg_list)
  output_final <- list(occ_partitions = occ_partitions, 
                       bg_partitions = bg_partitions)
  return(output_final)
}

# maxent_3D function outputs a list of models given what feature classes, 
# regularization multipliers, and partitions are tried.
# df_for_maxent = [[3]] output of brick_extract. dataframe where first column is vector
# of presences, other columns are environmental variables
# wanted_fc = feature classes to be tried in format c("L", "Q", "LQ", etc..)
# wanted_rm = regularization mutlipliers to be tried in format c(1:4)
# wanted_partition = output of partition_3D
# projection layers = list of spatraster stacks by depth for predictions
# occ_list = list of spatial points by depth of occurrences
maxent_3D <- function(df_for_maxent, wanted_fc, wanted_rm, wanted_partition,
                      projection_layers, occ_list) {
  # lets remove incomplete cases from the maxent_df so all points match
  maxent_df_present <- df_for_maxent %>% filter(presence == 1)
  maxent_df_absent <- df_for_maxent %>% filter(presence == 0)
  maxent_df_present <- maxent_df_present[complete.cases(maxent_df_present),]
  df_for_maxent <- rbind(maxent_df_present, maxent_df_absent)
 
  # first we'll make the fc input readable to maxent
  new_wanted_fc <- vector("list", length = length(wanted_fc))
  for(i in 1:length(wanted_fc)) {
    if(wanted_fc[i] == "L") {
        wanted_fc1 <- c("noquadratic", "noproduct", "nohinge", "noautofeature")
    } else if(wanted_fc[i] == "Q") {
        wanted_fc1 <- c("nolinear", "noproduct", "nohinge", "noautofeature")
    } else if(wanted_fc[i] == "H") {
        wanted_fc1 <- c("nolinear", "noquadratic", "noproduct", "noautofeature")
    } else if(wanted_fc[i] == "P") {
        wanted_fc1 <- c("nolinear", "noquadratic", "nohinge", "noautofeature")
    } else if(wanted_fc[i] == "LQ") {
        wanted_fc1 <- c("nohinge", "noproduct", "noautofeature")
    } else if(wanted_fc[i] == "LH") {
        wanted_fc1 <- c("noproduct", "noquadratic", "noautofeature")
    } else if(wanted_fc[i] == "LP") {
        wanted_fc1 <- c("noquadratic", "nohinge", "noautofeature")
    } else if(wanted_fc[i] == "QH") {
        wanted_fc1 <- c("nolinear", "noproduct", "noautofeature")
    } else if(wanted_fc[i] == "QP") {
        wanted_fc1 <- c("nolinear", "nohinge", "noautofeature")
    } else if(wanted_fc[i] == "HP") {
        wanted_fc1 <- c("nolinear", "noquadratic", "noautofeature")
    } else if(wanted_fc[i] == "LQP") {
        wanted_fc1 <- c("nohinge", "noautofeature")
    } else if(wanted_fc[i] == "LQH") {
        wanted_fc1 <- c("noproduct", "noautofeature")
    } else if(wanted_fc[i] == "QHP") {
        wanted_fc1 <- c("nolinear", "noautofeature")
    } else if(wanted_fc[i] == "LHP") {
        wanted_fc1 <- c("noquadratic", "noautofeature")
    } else if(wanted_fc[i] == "LQHP") {
        wanted_fc1 <- c("noautofeature")
    } else {
        print("not an applicable feature class")
    }
    new_wanted_fc[[i]] <- wanted_fc1
  }
  
  # next we'll make a full list of the fc and rm combinations
  # j loop for each rm, make a list the length of rm attached to fc[i]
  # so there will be i lists for each fc
  # concatenate the i lists into one list to make the full combined rm and fc list
  fc_perm_list <- vector("list", length = length(new_wanted_fc))
  for(i in 1:length(new_wanted_fc)) {
    rm_perm_list <- vector("list", length = length(wanted_rm)) 
      for(j in 1:length(wanted_rm)) {
        rm_element <- wanted_rm[j]
        fc_element <- new_wanted_fc[[i]]
        rm_perm_list[[j]] <- list(rm_element, fc_element)
      }
    fc_perm_list[[i]] <- rm_perm_list
  }
  final_perm_call <- vector("list", length = length(wanted_fc))
  for(i in 1:length(wanted_fc)) {
    rm_parts <- vector(length = length(wanted_rm))
    for(j in 1:length(wanted_rm)) {
      rm_parts[j] <- paste0("fc_perm_list[[", i, "]][[", j, "]]")
    }
    final_perm_call[[i]] <- rm_parts
  }
  final_perm_call_real <- unlist(final_perm_call)
  final_perm_list <- vector("list", length = length(final_perm_call_real))
  for(i in 1:length(final_perm_call_real)) {
    final_perm_list[[i]] <- eval(parse(text = final_perm_call_real[i]))
  }
  
  # if there is no partitioning scheme, a model will be run that is trained 
  # on all of the data. these models will also be returned as the models
  # run for a partitioning scheme, but the partitioning scheme run will also
  # return validation statistics. it can only return AUC
  if(is.null(wanted_partition)) {
    model_list <- vector("list", length = length(final_perm_list))
    auc_list <- vector("list", length = length(final_perm_list))
    total_coef_list <- vector(length = length(final_perm_list))
    nonzero_coef_list <- vector(length = length(final_perm_list))
    predictions_list <- vector("list", length = length(final_perm_list))
    AICc_list <- vector(length = length(final_perm_list))
    print("running models")
    for(i in 1:length(final_perm_list)) {
      # create model
      mod1 <- maxent(x = df_for_maxent[,-1], p = df_for_maxent[,1], 
                   args = c(paste0("betamultiplier=", 
                                   final_perm_list[[i]][[1]]),
                                   paste0(final_perm_list[[i]][[2]])))
      # save model
      model_list[[i]] <- mod1
      # retrieve AUC
      auc_list[[i]] <- mod1@results[which(rownames(mod1@results) == 
                                            "Training.AUC"),1]
      # retrieve coefficents
      all_coef <- as.numeric(unlist(strsplit(mod1@lambdas, ",")))
      all_coef <- all_coef[-which(is.na(all_coef))]
      total_coef_list[i] <- length(all_coef)
      nonzero_coef <- all_coef[which(all_coef != 0)]
      nonzero_coef <- nonzero_coef[which(nonzero_coef != 0.000000e+00)]
      nonzero_coef_list[i] <- length(nonzero_coef)
      # generate predictions and pull out values at occurrences for AICc calc
      predicted_suit_list <- vector("list", length = length(projection_layers))
      wanted_val_list <- vector("list", length = length(projection_layers))
      for(j in 1:length(projection_layers)) {
        predicted_suit <- predict(mod1, projection_layers[[j]])
        predicted_suit_list[[j]] <- predicted_suit
        standard_suit <- predicted_suit
        values(standard_suit) <- values(standard_suit)/max(values(standard_suit), 
                                                           na.rm = T)
        if(any(!(is.na(occ_list[[j]])))) {
          wanted_val_list[[j]] <- terra::extract(standard_suit, 
                                                 data.matrix(occ_list[[j]]))
        } else {
          wanted_val_list[[j]] <- NA
        }
      }
      wanted_val_final <- do.call(rbind, wanted_val_list)[,1]
      wanted_val_final <- wanted_val_final[which(!(is.na(wanted_val_final)))]
      predictions_list[[i]] <- rast(predicted_suit_list)
      # calculating AICc
      bigk <- length(nonzero_coef)
      likelihood <- sum(log(wanted_val_final))
      littlen <- length(wanted_val_final)
      if(bigk == (littlen - 1)) {
        AICc_list[i] <- NA
      } else {
        AICc <- ((2*bigk) - (2*likelihood)) + 
          (((2*bigk)*(bigk + 1))/(littlen - bigk - 1))
        AICc_list[i] <- AICc
      }
    }
    all_auc <- do.call(rbind, auc_list)
    rm <- rep(wanted_rm, times = length(wanted_fc))
    fc_list <- vector("list", length = length(wanted_fc))
    for(i in 1:length(wanted_fc)) {
      fc_list[[i]] <- rep(wanted_fc[i], times = length(wanted_rm))
    }
    fc <- unlist(fc_list)
    delta.AICc <- AICc_list - min(AICc_list, na.rm = T)
    mod_results <- data.frame(fc, rm, all_auc[,1], total_coef_list, nonzero_coef_list,
                              AICc_list, delta.AICc)
    colnames(mod_results) <- c("fc", "rm", "train.AUC", "toal.coef", "nonzero.coef",
                               "AICc", "delta.AICc")
    final_output <- list(models = model_list, predictions = predictions_list,
                         results = mod_results)
    
  # generating models for a given partition scheme  
  } else {
    
    # full models for model part of output
    model_list <- vector("list", length = length(final_perm_list))
    print("running models")
    for(i in 1:length(final_perm_list)) {
      # create model
      mod1 <- maxent(x = df_for_maxent[,-1], p = df_for_maxent[,1], 
                     args = c(paste0("betamultiplier=", 
                                     final_perm_list[[i]][[1]]),
                              paste0(final_perm_list[[i]][[2]])))
      # save model
      model_list[[i]] <- mod1
    }
    
    num_of_partitions <- unique(wanted_partition$bg_partitions)
    # generating training and testing models for each partition
    part_model_list <- vector("list", length = length(num_of_partitions))
    part_eval_list <- vector("list", length = length(num_of_partitions))
    part_predictions_list <- vector("list", length = length(num_of_partitions))
    present_only <- df_for_maxent[which(df_for_maxent$presence == 1),]
    absent_only <- df_for_maxent[which(df_for_maxent$presence == 0),]
    for (i in (1:length(num_of_partitions))) {
      test_p <- present_only[which(wanted_partition$occ_partitions == 
                                     num_of_partitions[i]),]
      train_p <- present_only[which(wanted_partition$occ_partitions != 
                                      num_of_partitions[i]),]
      test_a <- absent_only[which(wanted_partition$bg_partitions == 
                                    num_of_partitions[i]),]
      train_a <- absent_only[which(wanted_partition$bg_partitions != 
                                     num_of_partitions[i]),]
      test_full <- rbind(test_p, test_a)
      train_full <- rbind(train_p, train_a)
      mod_per_param_list <- vector("list", length = length(final_perm_list))
      auc_val_per_param_list <- vector(length = length(final_perm_list))
      auc_train_per_param_list <- vector("list", length = length(final_perm_list))
      total_coef_list <- vector(length = length(final_perm_list))
      nonzero_coef_list <- vector(length = length(final_perm_list))
      predictions_list <- vector("list", length = length(final_perm_list)) 
      print(paste0("training and testing models for partition ", i))
      for(j in 1:length(final_perm_list)) { 
        # generating each model combination for partition i
        train_model <- maxent(x = train_full[,-1], p = train_full[,1], 
                       args = c(paste0("betamultiplier=", 
                                       final_perm_list[[j]][[1]]),
                                paste0(final_perm_list[[j]][[2]])))
        # saving model
        mod_per_param_list[[j]] <- train_model
        # retrieving validation auc
        test_full_present <- test_full %>% filter(presence == 1)
        test_full_absent <- test_full %>% filter(presence == 0)
        ev <- dismo::evaluate(p = as.matrix(test_full_present[,-1]), 
                       a = as.matrix(test_full_absent[,-1]), 
                       model = train_model)
        auc_val_per_param_list[j] <- ev@auc
        # retrieving training AUC
        auc_train_per_param_list[[j]] <- 
          train_model@results[which(rownames(train_model@results) == "Training.AUC"),1]
        # retrieve coefficents
        all_coef <- as.numeric(unlist(strsplit(train_model@lambdas, ",")))
        all_coef <- all_coef[-which(is.na(all_coef))]
        total_coef_list[j] <- length(all_coef)
        nonzero_coef <- all_coef[which(all_coef != 0)]
        nonzero_coef <- nonzero_coef[which(nonzero_coef != 0.000000e+00)]
        nonzero_coef_list[j] <- length(nonzero_coef)
        # generate predictions; AICc will be calulated after model averaging across
        # partitions
        predicted_suit_list <- vector("list", length = length(projection_layers))
        wanted_val_list <- vector("list", length = length(projection_layers))
        for(k in 1:length(projection_layers)) {
          predicted_suit <- predict(train_model, projection_layers[[k]])
          predicted_suit_list[[k]] <- predicted_suit
        }
        predictions_list[[j]] <- predicted_suit_list
      }
      part_model_list[[i]] <- mod_per_param_list
      all_auc_train <- do.call(rbind, auc_train_per_param_list)
      rm <- rep(wanted_rm, times = length(wanted_fc))
      fc_list <- vector("list", length = length(wanted_fc))
      for(j in 1:length(wanted_fc)) {
        fc_list[[j]] <- rep(wanted_fc[j], times = length(wanted_rm))
      }
      fc <- unlist(fc_list)
      mod_results <- data.frame(fc, rm, all_auc_train, auc_val_per_param_list,
                              total_coef_list, 
                              nonzero_coef_list)
      colnames(mod_results) <- c("fc", "rm", "train.AUC", "val.AUC", "total.coef", 
                               "nonzero.coef")
      part_eval_list[[i]] <- mod_results
      part_predictions_list[[i]] <- predictions_list
    }
    
    # creating a mod results dataframe averaged across all partitions
    wantdf <- part_eval_list[[1]][,-c(1:2)]
    cell_num <- (nrow(wantdf)*ncol(wantdf))
    avg_mod_result <- matrix(nrow = nrow(wantdf), 
                             ncol = ncol(wantdf))
    for(i in 1:cell_num) {
      avg_vec <- vector(length = length(num_of_partitions))
      for(j in 1:length(num_of_partitions)) {
        wantdf <- part_eval_list[[j]][,-c(1:2)]
        avg_vec[j] <- as.matrix(wantdf)[i]
      }
      avg_mod_result[i] <- mean(avg_vec)
    }
    avg_mod_result <- data.frame(avg_mod_result)
    colnames(avg_mod_result) <- c("avg.train.AUC", "avg.val.AUC", "avg.total.coef", 
                                  "avg.nonzero.coef")
    avg_mod_result <- cbind(part_eval_list[[1]][,1:2], avg_mod_result)
    
    # creating an average suitability list across partitions for final predictions
    # and average AICc
    avg_predictions_list <- vector("list", length = length(final_perm_list))
    for(i in 1:length(final_perm_list)) {
      avg_suit_list <- vector("list", length = length(projection_layers))
      for(j in 1:length(projection_layers)) {
        partition_layers <- vector("list", length = length(num_of_partitions))
        for(k in 1:length(num_of_partitions)) {
          partition_layers[[k]] <- part_predictions_list[[k]][[i]][[j]]
        }
        avg_across_depth_for_model_i <- terra::app(rast(partition_layers), mean)
        avg_suit_list[[j]] <- avg_across_depth_for_model_i
      }
      avg_predictions_list[[i]] <- rast(avg_suit_list)
    }
    
    # calculating average AICc for each model type
    avg.AICc <- vector(length = length(final_perm_list))
    for(i in 1:length(final_perm_list)) {
      wanted_val_list <- vector("list", length = length(projection_layers))
      for(j in 1:length(projection_layers)) {
        standard_suit <- avg_predictions_list[[i]][[j]]
        values(standard_suit) <- values(standard_suit)/max(values(standard_suit), 
                                                         na.rm = T)
        if(any(!(is.na(occ_list[[j]])))) {
        wanted_val_list[[j]] <- terra::extract(standard_suit, 
                                               data.matrix(occ_list[[j]]))
        } else {
        wanted_val_list[[j]] <- NA
        }
      }
    wanted_val_final <- do.call(rbind, wanted_val_list)[,1]
    wanted_val_final <- wanted_val_final[which(!(is.na(wanted_val_final)))]
    # calculating AICc
    bigk <- avg_mod_result$avg.nonzero.coef[i]
    likelihood <- sum(log(wanted_val_final))
    littlen <- length(wanted_val_final)
    if(bigk == (littlen - 1)) {
      avg.AICc[i] <- NA
    } else {
      AICc <- ((2*bigk) - (2*likelihood)) + 
        (((2*bigk)*(bigk + 1))/(littlen - bigk - 1))
      avg.AICc[i] <- AICc
    }
    }
    avg.delta.AICc <- avg.AICc-min(avg.AICc)
    avg_mod_result <- cbind(avg_mod_result, avg.AICc, avg.delta.AICc)
    
    final_output <- list(models = model_list, predictions = avg_predictions_list,
                         results = avg_mod_result, partition_results = part_eval_list)
  }
  return(final_output)
}

# threshold_3D outputs a SpatRaster stack of thresholded rasters. It can either take
# the $predictions output of maxent_3D, which is the model projected back onto
# the accessible area, or a new projection layer set for the model to be projected
# and thresholded on.
# wanted_model = model for projecting
# projection_layers = list of rasterstacks by depth for the model to be projected onto
# predicted_layers = $predictions output of maxent_3D
# thresholding_vals = a vector of threshold percentile values
# occ_list = the list of occurrences by depth slice
# bg_list = the list of bachground values by depth slice
threshold_3D <- function(wanted_model, projection_layers, predicted_layers,
                         thresholding_vals, occ_list, bg_list) {
  # first for the case where the predictions from maxent_3D are provided
  if(is.na(projection_layers)) {
    sdm_suit_vals <- vector("list", length = length(occ_list))
    for(i in 1:length(occ_list)) {
      if(any(is.na(occ_list[[i]]))) {
        sdm_suit_vals[[i]] <- NA
      } else {
        sdm_suit_vals[[i]] <- terra::extract(predicted_layers[[i]], 
                        occ_list[[i]])
      }
    }
  sdm_suits <- do.call(rbind, sdm_suit_vals)[,2]
  sdm_suits <- sdm_suits[which(!(is.na(sdm_suits)))]
  sdm_suits <- sdm_suits[order(sdm_suits, decreasing = T)]
  thresh_sdm_list <- vector("list", length = length(thresholding_vals))
  tss_list <- vector(length = length(thresholding_vals))
  for(i in seq_along(thresholding_vals)) {
    thresholded_layers <- vector("list", length = length(predicted_layers))
    real_abs_list <- vector("list", length = length(bg_list))
    thresh <- round(length(sdm_suits)*thresholding_vals[i])
    thresh_val <- sdm_suits[thresh]
    for(j in 1:dim(predicted_layers)[3]) {
      thresholded_layers[[j]] <- rast(reclassify(raster(predicted_layers[[j]]), 
                                          rcl = c(0, thresh_val, 0, 
                                                    thresh_val, 1, 1)))
      real_abs_list[[j]] <- terra::extract(thresholded_layers[[j]], bg_list[[j]])
    }
    thresh_sdm_list[[i]] <- rast(thresholded_layers)
    real_abs <- do.call(rbind, real_abs_list)[,2]
    specificity <- length(which(real_abs == 0))/length(real_abs)
    sensitivity <- thresholding_vals[i]
    tss <- (sensitivity + ((1/3)*specificity)) - 1
    tss_list[i] <- tss
  }
  choose_final <- which(tss_list == min(tss_list))
  tss_df <- data.frame(thresholding_vals, tss_list)
  colnames(tss_df) <- c("threshold_value", "tss")
  final_threshold_output <- list(threshold_layers = thresh_sdm_list[[choose_final]],
                                 tss_results = tss_df)
  }
  return(final_threshold_output)
}

