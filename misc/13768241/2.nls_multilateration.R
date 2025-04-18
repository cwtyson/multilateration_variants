## CALIBRATION OF CTT DATA ##
# Multilateration using Nearest Least Squares in R

library(tidyverse)
library(sp)
library(sf)
library(lubridate)
library(cowplot)

## 1. Calculate estimated distances from decay model -----------------------------------

## AVG rss
# load dataframe with calculated distances, nodes, avg rss values per reloc
df_relocs <- read_csv("distance_curve_filtered_4.33.csv")

# calculated estimated distance from fitted general nls model
# use exponential model formula: avgRSS ~ a * exp(-S * distance) + K
# a = intercept
# S = decay factor
# K = horizontal asymptote

a <- 27.54383
S <- 0.008444905
K <- (-102.9043)

df_relocs$e_dist <- (log(df_relocs$avgRSS - K) - log(a)) / -S

## MAX rss
# load dataframe with calculated distances, nodes, max rss values per reloc
df_relocs <- read_csv("distance_curve_raw.csv")

# calculated estimated distance from fitted general nls model
# use exponential model formula: avgRSS ~ a * exp(-S * distance) + K
# a = intercept
# S = decay factor
# K = horizontal asymptote
# calculated from full model in 2.decay_function_time_diff
# final parameters for maxRSS model: a = 28.93531, S = 0.008566478, K = -102.4544

a <- 28.93531
S <- 0.008566478
K <- (-102.4544)

df_relocs$e_dist <- (log(df_relocs$maxRSS - K) - log(a)) / -S


## 2. Filter relocs with low signal or few nodes -----------------------------------

# remove NAs from df, when signal was lower than K
df_relocs <- df_relocs %>% filter(!is.na(e_dist))

# negative distances are assumed to be 0 (too close to node)
df_relocs <- df_relocs %>% mutate(e_dist = case_when(e_dist < 0 ~ 0,
                                                     e_dist >=0 ~ e_dist))

# cannot use relocations with less than 3 nodes
n_nodes_relocs <- df_relocs %>% group_by(reloc_id) %>% 
  summarise(n_nodes = n()) %>%
  ungroup()

relocs_over3_nodes <- n_nodes_relocs[n_nodes_relocs$n_nodes >2,]$reloc_id

df_relocs <- df_relocs %>% filter(reloc_id %in% relocs_over3_nodes)

# create vector of unique relocations
unique_relocations <- unique(df_relocs$reloc_id) # (over K and with over 3 nodes, 1772)


## 3. NLS multilateration: function ------------------------------------------------------

# Multilateration function, adapted from Paxton 2022

## obs! in all multilat functions, change avgRSS for maxRSS if
# running for maximum values

multilat <- do.call("rbind", lapply(unique_relocations, function(reloc) {
  
  # Separate by relocations
  reloc_single <- df_relocs %>% filter(reloc_id == reloc)
  
  # Calculate no nodes for the test
  no.nodes <- dplyr::n_distinct(reloc_single$NodeId)
  
  # Determine the node with the strongest RSS value
  max.RSS <- reloc_single %>% slice_max(avgRSS, n = 1, with_ties = FALSE)
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(reloc_id = character(), No.Nodes = numeric(), x.est = numeric(), y.est = numeric())
  
  # Optimize the estimate using least squares
      nls.test <- tryCatch(
    
    nls(e_dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = F, allpairs = T),
        data = reloc_single, start=list(NodeUTMx_solution=max.RSS$NodeUTMx, NodeUTMy_solution=max.RSS$NodeUTMy),
        control=nls.control(warnOnly = T, minFactor=1/1000, maxiter = 500)) ,
    error = function(e) NULL  # return NULL on error
    )
  
  if (!is.null(nls.test)) {
    # Determine an error around the point location estimate
    par.est <- cbind(coef(nls.test))
    
    # Estimated location of the point and error
    estimated.loc <- data.frame(
      reloc_id = reloc,
      No.Nodes = no.nodes,
      x.est = par.est[1, 1],
      y.est = par.est[2, 1]
    )
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
  }
  
  # Return dataframe
  estimated.location_results
  
}))


# join to multilat info
df_multilat_join <- df_relocs %>% select(reloc_id, lat, lon, TagId, elevation, tag_type, test_height, TestUTMx, TestUTMy, unique_trial_id) %>%
  distinct()

multilat_results <- left_join(multilat, df_multilat_join, by = "reloc_id")

# calculate difference between known and estimated locations
multilat_results$error_m <- raster::pointDistance(multilat_results[,c("x.est", "y.est")], multilat_results[,c("TestUTMx", "TestUTMy")], lonlat = F, allpairs = F)

## 4. NLS trilateration 3 strongest RSSI nodes: function -----------------------------------

trilat_rssi <- do.call("rbind", lapply(unique_relocations, function(reloc) {
  
  # Separate by relocations
  reloc_single <- df_relocs %>% filter(reloc_id == reloc)
  
  # Get only 3 top nodes
  reloc_top_rssi <- reloc_single %>% slice_max(avgRSS, n = 3, with_ties = FALSE) # only 3, no ties
  
  # Get location of nodes with signal for relocation
  node_UTM <- reloc_top_rssi %>% select(NodeUTMx, NodeUTMy) %>% distinct()
  
  # Calculate no nodes for the test
  no.nodes <- dplyr::n_distinct(reloc_top_rssi$NodeId)
  
  # Determine the node with the strongest RSS value
  max.RSS <- reloc_top_rssi %>% slice_max(avgRSS, n = 1, with_ties = FALSE)
  
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(reloc_id = character(), No.Nodes = numeric(), x.est = numeric(), y.est = numeric())
  
  
  # Optimize the estimate using least squares
  nls.test <- tryCatch(
    
    nls(e_dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = F, allpairs = T),
        data = reloc_top_rssi, start=list(NodeUTMx_solution=max.RSS$NodeUTMx, NodeUTMy_solution=max.RSS$NodeUTMy),
        control=nls.control(warnOnly = T, minFactor=1/1000, maxiter = 500)) ,
    error = function(e) NULL  # return NULL on error
  )
  
  if (!is.null(nls.test)) {
    # Determine an error around the point location estimate
    par.est <- cbind(coef(nls.test))
    
    # Estimated location of the point and error
    estimated.loc <- data.frame(
      reloc_id = reloc,
      No.Nodes = no.nodes,
      x.est = par.est[1, 1],
      y.est = par.est[2, 1]
    )
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
  }
  
  # Return dataframe
  estimated.location_results
  
}))

# join info
df_trilat_join <- df_relocs %>% select(reloc_id, lat, lon, TagId, elevation, tag_type, test_height, TestUTMx, TestUTMy, unique_trial_id) %>%
  distinct()

trilat_rssi_results <- left_join(trilat_rssi, df_trilat_join, by = "reloc_id")

# calculate difference between known and estimated locations
trilat_rssi_results$error_m <- raster::pointDistance(trilat_rssi_results[,c("x.est", "y.est")], trilat_rssi_results[,c("TestUTMx", "TestUTMy")], lonlat = F, allpairs = F)

## 5. NLS trilateration 3 closest nodes: function -----------------------------------

trilat_dist <- do.call("rbind", lapply(unique_relocations, function(reloc) {
  
  # Separate by relocations
  reloc_single <- df_relocs %>% filter(reloc_id == reloc)
  
  # Get only 3 top nodes
  reloc_top_dist <- reloc_single %>% slice_max(e_dist, n = 3, with_ties = FALSE) # only 3, no ties
  
  # Get location of nodes with signal for relocation
  node_UTM <- reloc_top_dist %>% select(NodeUTMx, NodeUTMy) %>% distinct()
  
  # Calculate no nodes for the test
  no.nodes <- dplyr::n_distinct(reloc_top_dist$NodeId)
  
  # Determine the node with the strongest RSS value
  max.RSS <- reloc_top_dist %>% slice_max(avgRSS, n = 1, with_ties = FALSE)
  
  
  # Create a dataframe for output estimates
  estimated.location_results <- data.frame(reloc_id = character(), No.Nodes = numeric(), x.est = numeric(), y.est = numeric())
  
  
  # Optimize the estimate using least squares
  nls.test <- tryCatch(
    
    nls(e_dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = F, allpairs = T),
        data = reloc_top_dist, start=list(NodeUTMx_solution=max.RSS$NodeUTMx, NodeUTMy_solution=max.RSS$NodeUTMy),
        control=nls.control(warnOnly = T, minFactor=1/1000, maxiter = 500)) ,
    error = function(e) NULL  # return NULL on error
  )
  
  if (!is.null(nls.test)) {
    # Determine an error around the point location estimate
    par.est <- cbind(coef(nls.test))
    
    # Estimated location of the point and error
    estimated.loc <- data.frame(
      reloc_id = reloc,
      No.Nodes = no.nodes,
      x.est = par.est[1, 1],
      y.est = par.est[2, 1]
    )
    
    # Populate dataframe with results
    estimated.location_results <- rbind(estimated.location_results, estimated.loc)
  }
  
  # Return dataframe
  estimated.location_results
  
}))

# join info
trilat_dist_results <- left_join(trilat_dist, df_trilat_join, by = "reloc_id")

# calculate difference between known and estimated locations
trilat_dist_results$error_m <- raster::pointDistance(trilat_dist_results[,c("x.est", "y.est")], trilat_dist_results[,c("TestUTMx", "TestUTMy")], lonlat = F, allpairs = F)

## 6. NLS all possible combinations for k-means clustering -----------------------------------

# multilat with clustering
multilat_kmeans <- do.call("rbind", lapply(unique_relocations, function(reloc) {
  
  cat("Working on iteration:", round((match(reloc, unique_relocations)/length(unique_relocations)*100), 2), "%\n")
  
  # Separate by relocations
  reloc_single <- df_trial %>% filter(reloc_id == reloc)
  
  # Maximum of 7 nodes with highest signal, for faster computing time
  reloc_single_sel <- reloc_single %>% slice_max(avgRSS, n = 7, with_ties = FALSE)
  
  # List nodes with detections in relocation
  nodes_reloc <- unique(reloc_single_sel$NodeId)
  
  # number of selected nodes
  no.nodes <- length(nodes_reloc)
  
  
  # Generate combinations for group sizes ranging from 3 to n
  all_combinations <- list()
  for (group_size in 3:no.nodes) {
    all_combinations[[as.character(group_size)]] <- combn(nodes_reloc, group_size)
  }
  
  # List of dataframes
  df_list <- names(all_combinations)
  
  # Perform trilateration for each combination
  estimated.location_results <- lapply(df_list, function(df) {
    
    df_sel <- as.data.frame(all_combinations[[df]])
    
    columns <- seq(1:ncol(df_sel))
    
    lapply(columns, function(c) {
      
      node_group <- df_sel[,c]
      
      reloc_3nodes <- reloc_single_sel %>% filter(NodeId %in% node_group)
      
      # Determine the node with the strongest RSS value
      max.RSS <- reloc_3nodes %>% slice_max(avgRSS, n = 1, with_ties = FALSE)
      
      # Optimize the estimate using least squares
      nls.test <- tryCatch(
        nls(e_dist ~ raster::pointDistance(data.frame(NodeUTMx, NodeUTMy), c(NodeUTMx_solution, NodeUTMy_solution), lonlat = FALSE, allpairs = TRUE),
            data = reloc_3nodes, start = list(NodeUTMx_solution = max.RSS$NodeUTMx, NodeUTMy_solution = max.RSS$NodeUTMy),
            control = nls.control(warnOnly = TRUE, minFactor = 1/1000, maxiter = 500)),
        error = function(e) NULL  # return NULL on error
      )
      
      if (!is.null(nls.test)) {
        # Get estimated x and y
        par.est <- cbind(coef(nls.test))
        
        # Estimated location of the point 
        data.frame(
          x.est = par.est[1, 1],
          y.est = par.est[2, 1],
          reloc_id = reloc_single$reloc_id,
          unique_trial_id = reloc_single$unique_trial_id,
          nodes_per_cluster = nrow(reloc_3nodes),
          original_no.nodes = length(unique(reloc_single$NodeId)),
          total_no.nodes = no.nodes,
          avg_cluster_accuracy = mean(reloc_3nodes$mean_accuracy)
        ) %>% distinct()
      } 
    })
  })
  
  # Transform to single dataframe
  loc_df <- bind_rows(map_dfr(estimated.location_results, bind_rows))
  
  
  # Return dataframe
  loc_df
  
}
))


## find optimal k within 95% ci
# within sum of squares reduction
# remove relocs with only 3 nodes detection = 1 possible loc estimate
only3nodes <- multilat_kmeans %>% group_by(reloc_id) %>% 
  summarise(no.est_relocs = n()) %>% 
  filter(no.est_relocs < 2)

only3nodes_relocs <- only3nodes$reloc_id

df_3nodes_relocs <- multilat_kmeans %>% filter(reloc_id %in% only3nodes_relocs)
df_sel <- multilat_kmeans %>% filter(!reloc_id %in% only3nodes_relocs)

# separate by relocations
unique_relocations <- unique(df_sel$reloc_id)

# wss function outside of loop
wss_function <- function(data, max_k = no.est_relocs -1) {
  wss <- sapply(1:max_k, function(k) {
    kmeans(data, centers = k)$tot.withinss
  })
  return(wss)
}


# function
optim_k_asymptote <- do.call("rbind", lapply(unique_relocations, function(reloc) {
  
  # separate by reloc
  reloc_single <- df_sel %>% filter(reloc_id == reloc)
  
  # Number of estimates
  no.est_relocs <- nrow(reloc_single)
  
  ### Find optimal k clusters
  kmeans_result <- NULL
  
  if (no.est_relocs > 5) {
    
    # run wss function
    wss_values <- wss_function(reloc_single[,c("x.est","y.est")], max_k = no.est_relocs -1)
    
    # set up dataframe
    wss_df <- data.frame(k = seq(1,length(wss_values)),
                         wss = wss_values)
    
    # Fit an exponential decay model to varying values of k
    decay_model <- tryCatch(nls(wss ~ a * exp(b * k) + z, 
                                start = list(a = max(wss_values), b = -0.1, z = min(wss_values)), 
                                data = wss_df),
                            error = function(e) NULL)  # return NULL on error
    
    if (!is.null(decay_model)) {
      # Extract the coefficient for the asymptotic value
      asymptote <- coef(decay_model)[3]
      
      # Define the percentage threshold (e.g., 95% of asymptotic value)
      ci <- tryCatch(confint(decay_model),
                     error = function(e) NULL)
      
      if (!is.null(ci)) {
        low_asym <- ci[3,1]
        hi_asym <- ci[3,2]
        
        # Select the lowest K within 95% of calculated asymptote
        wss_sel <- wss_df %>% filter(wss >= low_asym & wss <= hi_asym)
        optimal_k <- wss_sel[which.min(wss_sel$k),1]
        
        if (!is.null(optimal_k) && length(optimal_k) == 1) {
          
          ## Run k means clustering with optimal k
          kmeans_result <- kmeans(reloc_single[,c("x.est","y.est")], 
                                  algorithm = "Hartigan-Wong",
                                  centers = optimal_k,
                                  nstart = 10)
          
        } else { # if optimal k not found within 95% asymptote
          optimal_k <- 3
          
          ## Run k means clustering with optimal k
          kmeans_result <- kmeans(reloc_single[,c("x.est","y.est")], 
                                  algorithm = "Hartigan-Wong",
                                  centers = 3,
                                  nstart = 10)
        }
        
      } else { # if confidence intervals were not calculated
        optimal_k <- 3
        
        ## Run k means clustering with optimal k
        kmeans_result <- kmeans(reloc_single[,c("x.est","y.est")], 
                                algorithm = "Hartigan-Wong",
                                centers = 3,
                                nstart = 10)
      }
      
    } else {
      NULL  # Return NULL if decay_model is NULL for no.est_relocs >5
    }
    
  } else {
    if (no.est_relocs == 5) {
      ## Run k means clustering for 3 centers as a rule of thumb
      optimal_k <- 3
      
      kmeans_result <- kmeans(reloc_single[,c("x.est","y.est")], 
                              algorithm = "Hartigan-Wong",
                              centers = 3,
                              nstart = 10)
    } 
  }
  
  # Check if kmeans_result is not NULL and exists
  if (!is.null(kmeans_result) && exists("kmeans_result",  envir = environment())) {
    
    # df
    result_df <- data.frame(reloc_id = unique(reloc_single$reloc_id),
                            unique_trial_id = unique(reloc_single$unique_trial_id),
                            total_no.nodes = unique(reloc_single$total_no.nodes),
                            total_no.loc_estimates = no.est_relocs,
                            max_k = case_when(no.est_relocs > 1 ~ no.est_relocs - 1,
                                              no.est_relocs == 1 ~ 1),
                            optimal_k = optimal_k,
                            cluster_id = rownames(kmeans_result[["centers"]]),
                            cluster_x = kmeans_result[["centers"]][,1],
                            cluster_y = kmeans_result[["centers"]][,2],
                            cluster_no.relocations = kmeans_result$size)
    
    result_df
  } else {
    NULL  # Return NULL if kmeans_result is NULL
  }
  
}))



