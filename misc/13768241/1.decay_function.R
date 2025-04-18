## CALIBRATION OF CTT DATA ##
# Calculate exponential decay relation between
# signal strength and distance

# load packages
library(tidyverse)
library(sp)
library(sf)
library(lubridate)
library(nlstools)
library(boot)
library(cowplot)


## 1. Load and tidy data ----------------------------------------------
# import data
data_df <- read.csv("time_diff_matched_df.csv")

# give time format
data_df$time_difference_c<- gsub("S","", data_df$time_difference)
data_df$time_difference_c <- seconds(data_df$time_difference_c)
data_df$time_difference <- data_df$time_difference_c

data_df$abs_time_difference_c<- gsub("S","", data_df$abs_time_difference)
data_df$abs_time_difference_c <- seconds(data_df$abs_time_difference_c)
data_df$abs_time_difference <- data_df$abs_time_difference_c

data_df <- data_df %>% select(-c(time_difference_c, abs_time_difference_c))

# filter time differences larger than 10s
data_filtered <- data_df %>% filter(abs_time_difference < 10)

# calculate avgRSSI and maxRSSI by reloc
avg_rssi <- data_filtered %>% group_by(NodeId, reloc_id) %>%
  summarise(maxRSS = max(TagRSSI),
            avgRSS = mean(TagRSSI),
            sdRSS = sd(TagRSSI),
            n.det = n()) %>%
  ungroup()

quantile(avg_rssi$sdRSS, probs = 0.95, na.rm = TRUE) # some sd nas if there was only one det by node per reloc
quantile(avg_rssi$sdRSS, na.rm = TRUE) # 75% is 3.6 sd
max(avg_rssi$sdRSS, na.rm = TRUE) # 21.2132

quantile(avg_rssi$n.det)

# give mean tagrssi per node and reloc
dat <- left_join(data_filtered, avg_rssi, by = c("NodeId", "reloc_id"))

unique(is.na(dat))

# calculate difference from mean of Node and reloc
dat$avgRSS_dif <- abs(dat$avgRSS - dat$TagRSSI)
quantile(dat$avgRSS_dif, probs = 0.95) # 4.33
quantile(dat$avgRSS_dif)
max(dat$avgRSS_dif) # 15.67

# remove any points further than 4.33 away from mean (95% sd)
dat2 <- dat %>% filter(avgRSS_dif < 4.33)
quantile(dat2$avgRSS_dif, probs = 0.95) # 3
max(dat2$avgRSS_dif) # 4.25


# calculate avg rssi with new filtered data
avg_rssi2 <- dat2 %>% group_by(NodeId, reloc_id) %>%
  summarise(avgRSS = mean(TagRSSI),
            sdRSS = sd(TagRSSI),
            n.det = n()) %>%
  ungroup()


quantile(avg_rssi2$sdRSS, probs = 0.95, na.rm = TRUE) # 4.95
range(avg_rssi2$sdRSS, na.rm = TRUE)

# join info 
join_df <- dat %>% dplyr::select(-c(time_difference, abs_time_difference, TagRSSI,
                                     Time.local, time_format, node_local_time,
                                     avgRSS_dif, maxRSS, avgRSS, sdRSS, n.det # remove also first averages calculated
                                     )) %>%
  distinct()

df2 <- left_join(avg_rssi2, join_df, by = c("NodeId", "reloc_id"))

unique(is.na(df2)) # check join

# for maxRSS
df <- dat %>% dplyr::select(-c(time_difference, abs_time_difference, TagRSSI,
                               Time.local, time_format, node_local_time,
                               avgRSS_dif, avgRSS))
df <- df %>% distinct()

## 2. Calculate distance between known locations and all nodes ------------------------------------

# import node data
nodes <- read.csv("Nodes.csv")

# df no filter signals for each NodeId and reloc_id
# df2 filtered with dif below 4.62

# spatial coordinates
original_crs <- st_crs("EPSG:4326") # usual lon lat
new_crs <- st_crs("EPSG:3116") # Bogota Zone

# convert test coordinates into UTM, projected
# df
sf_df <- st_as_sf(df, coords = c("lon", "lat"), crs = original_crs)
sf_df_transformed <- st_transform(sf_df, crs = new_crs)
coordinates <- st_coordinates(sf_df_transformed)

# convert lon into Easting location
df$TestUTMx <- coordinates[,"X"]

# convert lat into Northing location
df$TestUTMy <- coordinates[,"Y"]

# create unique row id
df$row_id <- seq(1:length(df$TagId))

# df2
# df
sf_df2 <- st_as_sf(df2, coords = c("lon", "lat"), crs = original_crs)
sf_df_transformed2 <- st_transform(sf_df2, crs = new_crs)
coordinates2 <- st_coordinates(sf_df_transformed2)

# convert lon into Easting location
df2$TestUTMx <- coordinates2[,"X"]

# convert lat into Northing location
df2$TestUTMy <- coordinates2[,"Y"]

# create unique row id
df2$row_id <- seq(1:length(df2$TagId))


## for nodes too, double checking previously recorded in 1.Nodes_files_prep
# df
sf_nodes <- st_as_sf(nodes, coords = c("Longitude", "Latitude"), crs = original_crs)
nodes_transformed <- st_transform(sf_nodes, crs = new_crs)
coordinates_nodes <- st_coordinates(nodes_transformed)

# convert lon into Easting location
nodes$TestUTMx_c <- coordinates_nodes[,"X"]

# convert lat into Northing location
nodes$TestUTMy_c <- coordinates_nodes[,"Y"]


# calculate all pairwise distances
dst1 <- raster::pointDistance(df[,c("TestUTMx", "TestUTMy")], nodes[,c("NodeUTMx", "NodeUTMy")], lonlat = F, allpairs = T)
dst2 <- raster::pointDistance(df2[,c("TestUTMx", "TestUTMy")], nodes[,c("NodeUTMx", "NodeUTMy")], lonlat = F, allpairs = T)

# Make matrix into a dataframe
dist_df1 <- data.frame(dst1, row.names = df$row_id)
colnames(dist_df1) <- nodes$NodeId
dist_df1$row_id <- rownames(dist_df1)

dist_df2 <- data.frame(dst2, row.names = df2$row_id)
colnames(dist_df2) <- nodes$NodeId
dist_df2$row_id <- rownames(dist_df2)

# rearrange data
dist_df_tidy1 <- dist_df1 %>%
  tidyr::gather(key = "NodeId", value = "distance", -row_id)
dist_df_tidy2 <- dist_df2 %>%
  tidyr::gather(key = "NodeId", value = "distance", -row_id)

# add relevant info
test_info <- df %>% dplyr::select(row_id, NodeId, TagId, maxRSS, sdRSS, n.det, elevation,
                                       trial_num, lon, lat, tag_type, test_local_time, start_time, end_time,
                                       test_height, tag_type, TestUTMx, TestUTMy, unique_trial_id,
                                       reloc_id)
test_info2 <- df2 %>% dplyr::select(row_id, NodeId, TagId, avgRSS, sdRSS, n.det, elevation,
                                  trial_num, lon, lat, tag_type, test_local_time, start_time, end_time,
                                  test_height, tag_type, TestUTMx, TestUTMy, unique_trial_id,
                                  reloc_id)

node_info <- nodes %>% dplyr::select(NodeId, Latitude, Longitude,
                                     NodeUTMx, NodeUTMy) %>%
  rename(node_lon = Longitude, node_lat = Latitude)

test_info$row_id <- as.character(test_info$row_id)
test_info2$row_id <- as.character(test_info2$row_id)

# CAREFUL! Node 3288e6 as 3288E6 in dist_df_tidy and node_info! 3798b0 not in filtered
base::setdiff(dist_df_tidy1$NodeId, test_info$NodeId)
base::setdiff(node_info$NodeId, test_info$NodeId) 
dist_df_tidy1[dist_df_tidy1$NodeId == "3288E6",]$NodeId <- "3288e6"
dist_df_tidy2[dist_df_tidy2$NodeId == "3288E6",]$NodeId <- "3288e6"

node_info[node_info$NodeId == "3288E6",]$NodeId <- "3288e6"

dist_df_full1 <- left_join(test_info, dist_df_tidy1, by = c("row_id", "NodeId"))
dist_df_full1 <- left_join(dist_df_full1, node_info, by = "NodeId")

dist_df_full2 <- left_join(test_info2, dist_df_tidy2, by = c("row_id", "NodeId"))
dist_df_full2 <- left_join(dist_df_full2, node_info, by = "NodeId")

unique(is.na(dist_df_full1))
unique(is.na(dist_df_full2))


# save
write_csv(dist_df_full1, "distance_curve_raw.csv") # df no filter signals for each NodeId and reloc_id
write_csv(dist_df_full2, "distance_curve_filtered_4.33.csv") # df2 filtered with dif below 4.33


## 3. Fit decay model with nls ------------------------------------

# read again
dist_df_full1 <- read.csv("distance_curve_raw.csv") # df no filter signals for each NodeId and reloc_id
dist_df_full2 <- read.csv("distance_curve_filtered_4.33.csv") # df2 filtered with dif below 4.33

df_model1 <- dist_df_full1
df_model2 <- dist_df_full2

# maxRSS, choosing max signal of same "beep": reloc_id, NodeId
# fit with horizontal distance only
# SSasymp is a asymptotic regression self-starting model
exp.mod1 <- nls(maxRSS ~ SSasymp(distance, Asym, R0, lrc), data = df_model1)
summary(exp.mod1)

exp(coef(exp.mod1)[["lrc"]])

# exponential model formula: RSS ~ a * exp(-S * distance) + K
# a = intercept
# S = decay factor
# K = horizontal asymptote
# values for general initial model, filtered data
a <- -73.51906
S <- 0.008566477 # exp(coef(exp.mod)[["lrc"]])
K <- -102.45436

nls.mod1 <- nls(maxRSS ~ a * exp(-S * distance) + K, start = list(a = a, S = S, K= K), 
                data = df_model1)
summary(nls.mod1)
coef(nls.mod1)
confint(nls.mod1)

# calculate boostrap CI
boot1 <- nlsBoot(nls.mod1)
summary(boot1)

# Get SE and save into df
se_mod1 <- summary(nls.mod1)$sigma
df_model1$se <- se_mod1

# final parameters for maxRSS model: a = 28.93531, S = 0.008566322, K = -102.4544

# Get residuals
df_model1$resid <- residuals(nls.mod1)

# Model rank
AIC(nls.mod1)
BIC(nls.mod1)
logLik(nls.mod1)

# Get model predictions
df_model1$pred <- predict(nls.mod1)

# Graph
(mod1_g <- ggplot() + 
    geom_point(data = df_model1, aes(x = distance, y = maxRSS), size = 1, alpha = 0.5, color = "#454545") +
    geom_ribbon(data = df_model1, aes(x = distance, ymin = pred - se, ymax = pred + se), fill = "#689FBB", alpha = 0.3) +
    geom_line(data = df_model1, aes(x = distance, y = pred), color = "#689FBB") +
    labs(x = "Distance (m)", y = "Maximum RSSI", color = "") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10)))

# avgRSS, filtered with 4.33 db max dif same "beep": reloc_id, NodeId
# fit with horizontal distance only
# SSasymp is a asymptotic regression self-starting model
exp.mod2 <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc), data = df_model2)

summary(exp.mod2)

exp(coef(exp.mod2)[["lrc"]])

# exponential model formula: avgRSS ~ a * exp(-S * distance) + K
# a = intercept
# S = decay factor
# K = horizontal asymptote
# values for general initial model, filtered data
a <- -75.36047
S <- 0.008444905 # exp(coef(exp.mod)[["lrc"]])
K <- -102.90431 

nls.mod2 <- nls(avgRSS ~ a * exp(-S * distance) + K, start = list(a = a, S = S, K= K), 
               data = df_model2)
summary(nls.mod2)
coef(nls.mod2)
confint(nls.mod2)

# calculate boostrap CI
boot2 <- nlsBoot(nls.mod2)
summary(boot2)

# Get SE and save into df
se_mod2 <- summary(nls.mod2)$sigma
df_model2$se <- se_mod2

# final parameters for filtered avgRSS model: a = 27.54339, S = 0.008444648, K = -102.9044

# Get residuals
df_model2$resid <- residuals(nls.mod2)

# Model rank
AIC(nls.mod2)
BIC(nls.mod2)
logLik(nls.mod2)

# Get model predictions
df_model2$pred <- predict(nls.mod2)

# Graph
(mod2_g <- ggplot() + 
    geom_point(data = df_model2, aes(x = distance, y = avgRSS), size = 1, alpha = 0.5, color = "#454545") +
    geom_ribbon(data = df_model2, aes(x = distance, ymin = pred - se, ymax = pred + se), fill = "#689FBB", alpha = 0.3) +
    geom_line(data = df_model2, aes(x = distance, y = pred), color = "#689FBB") +
    labs(x = "Distance (m)", y = "Average RSSI", color = "") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10)))


## 4. Fit decay models by node ------------------------------------
# repeat for df_model max and avg RSS
node_mod_list <- list()

node_ids <- unique(df_model$NodeId)
for (node_id in node_ids) {
  
  df_node <- df_model %>% filter(NodeId == node_id)
  
  if (nrow(df_node) > 0) {
    # Fit the model
    nls_model <- nls(avgRSS ~ a * exp(-S * distance) + K, 
                     start = list(a = a, S = S, K = K), 
                     data = df_node)
    
    # store results
    node_mod_results <- data.frame(NodeId = node_id,
                              n = nrow(df_node),
                              a = coef(nls_model)[1],
                              S = coef(nls_model)[2],
                              K = coef(nls_model)[3])
    
    node_mod_list[[node_id]] <- node_mod_results
  } else {
    # If there are not enough observations, store default values
    node_mod_results <- data.frame(NodeId = node_id,
                              n = 0,
                              a = a,
                              S = S,
                              K = K)
    node_mod_list[[node_id]] <- mod_results
  }
}

# Combine all model parameters
node_mod_df <- do.call(rbind, node_mod_list)

## 5. Fit decay models by tag ----------------------------------------------

# fit exp models, looping over each one
tags_sel <- unique(df_sel1$TagId)
base::setdiff(tags, tags_sel) # check removed
tags_sel_check <- unique(df_sel2$TagId) # check, must be the same


for (i in 1:length(tags_sel)) {
  
  tag <- tags_sel[i]
  
  ## for maxRSS
  # filter by tag
  df_1 <- df_sel1 %>% filter(TagId == tag)
  
  # fit initial, self-starting parameters
  exp.mod1 <- nls(maxRSS ~ SSasymp(distance, Asym, R0, lrc),
                  data = df_1)
  
  
  # exponential model formula: RSS ~ a * exp(-S * distance) + K
  # a = intercept
  # S = decay factor
  # K = horizontal asymptote
  a_maxRSS <- coef(exp.mod1)[["R0"]]
  S_maxRSS <- exp(coef(exp.mod1)[["lrc"]])
  K_maxRSS <- coef(exp.mod1)[["Asym"]]
  
  nls.mod1 <- nls(maxRSS ~ a_maxRSS * exp(-S_maxRSS * distance) + K_maxRSS,
                  start = list(a_maxRSS = a_maxRSS, S_maxRSS = S_maxRSS, K_maxRSS= K_maxRSS), 
                  data = df_1)
  
  
  # calculate boostrap CI
  boot1 <- nlsBoot(nls.mod1)
  
  # append to df
  df_1$a <- coef(nls.mod1)[["a_maxRSS"]]
  df_1$S <- exp(coef(nls.mod1)[["S_maxRSS"]])
  df_1$K <- coef(nls.mod1)[["K_maxRSS"]]
  df_1$se <- summary(nls.mod1)$sigma
  
  df_1$a_se <- boot1[["estiboot"]][1,2]
  df_1$a_2.5 <- boot1[["bootCI"]][1,2]
  df_1$a_97.5 <- boot1[["bootCI"]][1,3]
  
  df_1$S_se <- boot1[["estiboot"]][2,2]
  df_1$S_2.5 <- boot1[["bootCI"]][2,2]
  df_1$S_97.5 <- boot1[["bootCI"]][2,3]
  
  df_1$K_se <- boot1[["estiboot"]][3,2]
  df_1$K_2.5 <- boot1[["bootCI"]][3,2]
  df_1$K_97.5 <- boot1[["bootCI"]][3,3]
  
  # Get residuals
  df_1$resid <- residuals(nls.mod1)
  
  # Get model predictions
  df_1$pred <- predict(nls.mod1)
  
  # save
  write_csv(df_1, paste0("models_by_tag/no_smooth_maxrss/", tag, ".csv"))
  
  
  ## for avgRSS
  # filter by tag
  df_2 <- df_sel2 %>% filter(TagId == tag)
  
  # fit initial, self-starting parameters
  exp.mod2 <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc),
                  data = df_2)
  
  
  # exponential model formula: RSS ~ a * exp(-S * distance) + K
  # a = intercept
  # S = decay factor
  # K = horizontal asymptote
  a_avgRSS <- coef(exp.mod2)[["R0"]]
  S_avgRSS <- exp(coef(exp.mod2)[["lrc"]])
  K_avgRSS <- coef(exp.mod2)[["Asym"]]
  
  nls.mod2 <- nls(avgRSS ~ a_avgRSS * exp(-S_avgRSS * distance) + K_avgRSS,
                  start = list(a_avgRSS = a_avgRSS, S_avgRSS = S_avgRSS, K_avgRSS= K_avgRSS), 
                  data = df_2)
  
  
  # calculate boostrap CI
  boot2 <- nlsBoot(nls.mod2)
  
  # append to df
  df_2$a <- coef(nls.mod2)[["a_avgRSS"]]
  df_2$S <- exp(coef(nls.mod2)[["S_avgRSS"]])
  df_2$K <- coef(nls.mod2)[["K_avgRSS"]]
  df_2$se <- summary(nls.mod2)$sigma
  
  df_2$a_se <- boot2[["estiboot"]][1,2]
  df_2$a_2.5 <- boot2[["bootCI"]][1,2]
  df_2$a_97.5 <- boot2[["bootCI"]][1,3]
  
  df_2$S_se <- boot2[["estiboot"]][2,2]
  df_2$S_2.5 <- boot2[["bootCI"]][2,2]
  df_2$S_97.5 <- boot2[["bootCI"]][2,3]
  
  df_2$K_se <- boot2[["estiboot"]][3,2]
  df_2$K_2.5 <- boot2[["bootCI"]][3,2]
  df_2$K_97.5 <- boot2[["bootCI"]][3,3]
  
  # Get residuals
  df_2$resid <- residuals(nls.mod2)
  
  # Get model predictions
  df_2$pred <- predict(nls.mod2)
  
  # save
  write_csv(df_2, paste0("models_by_tag/no_smooth_avgrss/", tag, ".csv"))
  
  
}

