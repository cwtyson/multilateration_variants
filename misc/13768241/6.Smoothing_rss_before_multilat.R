## CALIBRATION OF CTT DATA ##

# Smoothing RSSI, as time series
# Kalman filters

library(tidyverse)
library(lubridate)
library(KFAS)
library(xts)


## 1. Import and tidy data -----------------------------------------------
## avgrssi
# import data
df_avgrssi <- read_csv("distance_curve_filtered_4.33.csv")

# check? time zone imported correctly?
time_check <- df_avgrssi

time_check$h <- hour(time_check$test_local_time)

(g_time <- ggplot() +
    geom_density(data = time_check, aes(x = h)))

# convert time zone!
df_avgrssi$test_local_time <- with_tz(df_avgrssi$test_local_time, 
                             tzone = "America/Bogota")

# check again
time_check <- df_avgrssi

time_check$h <- hour(time_check$test_local_time)

(g_time <- ggplot() +
    geom_density(data = time_check, aes(x = h)))


# calculate track duration for each point within track
# results in secs
df_avgrssi$duration <- df_avgrssi$test_local_time - df_avgrssi$start_time

# deviation of rssi from median
med_avgrssi <- df_avgrssi %>% group_by(unique_trial_id, NodeId) %>%
  summarise(med_avgrss = median(avgRSS))

j_avgrssi <- left_join(df_avgrssi, med_avgrssi, by = c("unique_trial_id", "NodeId"))

j_avgrssi <- j_avgrssi %>% mutate(dev_avgRSS = (avgRSS - med_avgrss) / med_avgrss)

# visualize
# all together
(dev_plot1 <- ggplot() +
    geom_line(data = j_avgrssi, aes(x = (duration/60), y = dev_avgRSS,
                                    group = interaction(unique_trial_id, NodeId),
                                    color = interaction(unique_trial_id, NodeId)),
              show.legend = FALSE, alpha = 0.3) +
    labs(y = "Deviation from median RSS", x = "Time elapsed (min)", title = "Without smoothing") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank()))


# separate by groups of tracks
unique_trials <- unique(j_avgrssi$unique_trial_id)

## max rssi
# import data
df_maxrssi <- read_csv("distance_curve_raw.csv")

# check? time zone imported correctly?
time_check <- df_maxrssi

time_check$h <- hour(time_check$test_local_time)

(g_time <- ggplot() +
    geom_density(data = time_check, aes(x = h)))

# convert time zone!
df_maxrssi$test_local_time <- with_tz(df_maxrssi$test_local_time, 
                                      tzone = "America/Bogota")

# check again
time_check <- df_maxrssi

time_check$h <- hour(time_check$test_local_time)

(g_time <- ggplot() +
    geom_density(data = time_check, aes(x = h)))


# calculate track duration for each point within track
# results in secs
df_maxrssi$duration <- df_maxrssi$test_local_time - df_maxrssi$start_time

# deviation of rssi from median
med_maxrssi <- df_maxrssi %>% group_by(unique_trial_id, NodeId) %>%
  summarise(med_maxrss = median(maxRSS))

j_maxrssi <- left_join(df_maxrssi, med_maxrssi, by = c("unique_trial_id", "NodeId"))

j_maxrssi <- j_maxrssi %>% mutate(dev_maxRSS = (maxRSS - med_maxrss) / med_maxrss)

# visualize
# all together
(dev_plot_max1 <- ggplot() +
    geom_line(data = j_maxrssi, aes(x = (duration/60), y = dev_maxRSS,
                                    group = interaction(unique_trial_id, NodeId),
                                    color = interaction(unique_trial_id, NodeId)),
              show.legend = FALSE, alpha = 0.3) +
    labs(y = "Deviation from median RSS", x = "Time elapsed (min)", title = "Without smoothing") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank()))


## 2. Kalman: avgRSS ---------------------------------------------------------
# run by node for each track
kalman_smoother_avgrss <- do.call("rbind", lapply(unique_trials, function(trial) {
  
track <- j_avgrssi %>% filter(unique_trial_id == trial)
nodes_in_track <- unique(track$NodeId)

ks_by_node <- do.call("rbind", lapply(nodes_in_track, function(node) {

node_seq <- track %>% filter(NodeId == node)

# make equal time intervals
ts <- xts(node_seq$avgRSS, order.by = node_seq$test_local_time,
          tzone = "America/Bogota")

# Make the time series regularly spaced with NAs for missing times
r_ts <- seq(start(ts), end(ts), by = "1 sec")
m_ts <- xts(merge(ts, xts(, order.by = r_ts)), order.by = r_ts)

# fit state space model, with non-informative initial parameters
ssm <- SSModel(m_ts$ts ~ SSMtrend(1, Q = list(matrix(NA), matrix(NA))), H = matrix(NA))
fit2 <- fitSSM(ssm, inits = c(0, 0, 0, 0), method = "BFGS")
ks_result <- KFS(fit2$model, filtering = "state", smoothing = "state")

# only get results for obs time stamps
non_na_indices <- which(!is.na(m_ts$ts))

k_filter_results <- ks_result$att
k_filter_sel <- k_filter_results[non_na_indices]

k_smooth_results <- ks_result$alphahat
k_smooth_sel <- k_smooth_results[non_na_indices]

data.frame(row_id = node_seq$row_id,
           NodeId = node_seq$NodeId,
           unique_trial_id = node_seq$unique_trial_id,
           reloc_id = node_seq$reloc_id,
           n.node_signals = length(node_seq$row_id),
           avgRSS_kf = as.vector(k_filter_sel), # filtered estimate of states
           avgRSS_ks = as.vector(k_smooth_sel)) # smoothed estimate of states
}))
}
))

# bind result to dataframe
kalman_results_avgrss <- left_join(j_avgrssi, kalman_smoother_avgrss, by = c("row_id", "NodeId",
                                                               "unique_trial_id", "reloc_id"))

# check join worked ok
unique(is.na(kalman_results_avgrss))

# save
write_csv(kalman_results_avgrss, "Kalman_smooth/tracks/avg_rss_kalman.csv")



## 3. Kalman: maxRSS ---------------------------------------------------------

# Apply Kalman Smoother
# separate by groups of tracks
unique_trials <- unique(j_maxrssi$unique_trial_id)

# run by node for each track
kalman_smoother_maxrss <- do.call("rbind", lapply(unique_trials, function(trial) {
  
  track <- j_maxrssi %>% filter(unique_trial_id == trial)
  nodes_in_track <- unique(track$NodeId)
  
  ks_by_node <- do.call("rbind", lapply(nodes_in_track, function(node) {
    
    node_seq <- track %>% filter(NodeId == node)
    
    # make equal time intervals
    ts <- xts(node_seq$maxRSS, order.by = node_seq$test_local_time,
              tzone = "America/Bogota")
    
    # Make the time series regularly spaced with NAs for missing times
    r_ts <- seq(start(ts), end(ts), by = "1 sec")
    m_ts <- xts(merge(ts, xts(, order.by = r_ts)), order.by = r_ts)
    
    
    # fit state space model, with non-informative initial parameters
    ssm <- SSModel(m_ts$ts ~ SSMtrend(1, Q = list(matrix(NA), matrix(NA))), H = matrix(NA))
    fit2 <- fitSSM(ssm, inits = c(0, 0, 0, 0), method = "BFGS")
    ks_result <- KFS(fit2$model, filtering = "state", smoothing = "state")
    
    # only get results for obs time stamps
    non_na_indices <- which(!is.na(m_ts$ts))
    
    k_filter_results <- ks_result$att
    k_filter_sel <- k_filter_results[non_na_indices]
    
    k_smooth_results <- ks_result$alphahat
    k_smooth_sel <- k_smooth_results[non_na_indices]
    
    data.frame(row_id = node_seq$row_id,
               NodeId = node_seq$NodeId,
               unique_trial_id = node_seq$unique_trial_id,
               reloc_id = node_seq$reloc_id,
               n.node_signals = length(node_seq$row_id),
               maxRSS_kf = as.vector(k_filter_sel), # filtered estimate of states
               maxRSS_ks = as.vector(k_smooth_sel)) # smoothed estimate of states
  }))
}
))

# bind result to dataframe
kalman_results_maxrss <- left_join(j_maxrssi, kalman_smoother_maxrss, by = c("row_id", "NodeId",
                                                                             "unique_trial_id", "reloc_id"))

# check join worked ok
unique(is.na(kalman_results_maxrss))

# save
write_csv(kalman_results_maxrss, "Kalman_smooth/tracks/max_rss_kalman.csv")

## 4. Cubic splines: avgRSS ---------------------------------------------------------

# run by node for each track
cubic_spline_avgrss <- do.call("rbind", lapply(unique_trials, function(trial) {
  
  track <- j_avgrssi %>% filter(unique_trial_id == trial)
  nodes_in_track <- unique(track$NodeId)
  
  ks_by_node <- do.call("rbind", lapply(nodes_in_track, function(node) {
    
    node_seq <- track %>% filter(NodeId == node)
    
    if (length(node_seq$row_id) >= 5) {
      
      sp <- smooth.spline(x = node_seq$test_local_time, y = node_seq$avgRSS)
      
      data.frame(NodeId = node_seq$NodeId,
                 unique_trial_id = node_seq$unique_trial_id,
                 reloc_id = node_seq$reloc_id,
                 n.node_signals = length(node_seq$row_id),
                 avgRSS_cs = sp[["y"]], # smoothed value
                 df_cs = sp[["df"]], # degrees of freedom
                 spar_cs = sp[["spar"]]) # smoothing determined by cross-validation
      
    } else {
      
      data.frame(NodeId = unique(node_seq$NodeId),
                 unique_trial_id = unique(node_seq$unique_trial_id),
                 reloc_id = unique(node_seq$reloc_id),
                 n.node_signals = length(node_seq$row_id),
                 avgRSS_cs = NA, 
                 df_cs = NA, 
                 spar_cs = NA)
      
    }
  }))
}))

# bind result to dataframe
cs_results_avgrss <- left_join(j_avgrssi, cubic_spline_avgrss, 
                               by = c("NodeId", "unique_trial_id", "reloc_id"))

# if less than 5 detections by node, use unsmoothed avgRSS
less_5 <- cs_results_avgrss %>% filter(is.na(avgRSS_cs))
range(less_5$n.node_signals)

cs_results_avgrss$avgRSS_cs <- ifelse(is.na(cs_results_avgrss$avgRSS_cs), cs_results_avgrss$avgRSS, cs_results_avgrss$avgRSS_cs)
unique(is.na(cs_results_avgrss))

# save
write_csv(cs_results_avgrss, "Cubic_spline_smooth/tracks/avg_rss_cs.csv")

## 5. Cubic splines: maxRSS ---------------------------------------------------------
# run by node for each track
cubic_spline_maxrss <- do.call("rbind", lapply(unique_trials, function(trial) {
  
  track <- j_maxrssi %>% filter(unique_trial_id == trial)
  nodes_in_track <- unique(track$NodeId)
  
  ks_by_node <- do.call("rbind", lapply(nodes_in_track, function(node) {
    
    node_seq <- track %>% filter(NodeId == node)
    
    if (length(node_seq$row_id) >= 5) {
      
      sp <- smooth.spline(x = node_seq$test_local_time, y = node_seq$maxRSS)
      
      data.frame(NodeId = node_seq$NodeId,
                 unique_trial_id = node_seq$unique_trial_id,
                 reloc_id = node_seq$reloc_id,
                 n.node_signals = length(node_seq$row_id),
                 maxRSS_cs = sp[["y"]], # smoothed value
                 df_cs = sp[["df"]], # degrees of freedom
                 spar_cs = sp[["spar"]]) # smoothing determined by cross-validation
      
    } else {
      
      data.frame(NodeId = unique(node_seq$NodeId),
                 unique_trial_id = unique(node_seq$unique_trial_id),
                 reloc_id = unique(node_seq$reloc_id),
                 n.node_signals = length(node_seq$row_id),
                 maxRSS_cs = NA, 
                 df_cs = NA, 
                 spar_cs = NA)
      
    }
  }))
}))

# bind result to dataframe
cs_results_maxrss <- left_join(j_maxrssi, cubic_spline_maxrss, 
                               by = c("NodeId", "unique_trial_id", "reloc_id"))

# if less than 5 detections by node, use unsmoothed maxRSS
less_5 <- cs_results_maxrss %>% filter(is.na(maxRSS_cs))
range(less_5$n.node_signals)

cs_results_maxrss$maxRSS_cs <- ifelse(is.na(cs_results_maxrss$maxRSS_cs), cs_results_maxrss$maxRSS, cs_results_maxrss$maxRSS_cs)
unique(is.na(cs_results_maxrss))

# save
write_csv(cs_results_maxrss, "Cubic_spline_smooth/tracks/max_rss_cs.csv")
