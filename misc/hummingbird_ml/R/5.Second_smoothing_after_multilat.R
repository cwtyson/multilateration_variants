## CALIBRATION OF CTT DATA ##
# Smooth location estimates

# load required packages
library(tidyverse)
library(KFAS)
library(xts)

## Import localization data -----------------------------------------------------

df_all <- read.csv("df_all_models_error.csv")

unique_loc <- unique(df_all$loc_method)
unique_param <- unique(df_all$param)
unique_signal <- unique(df_all$type_signal)
unique_smooth <- unique(df_all$smoothing)

method_comb <- expand.grid(
  loc_method = unique_loc,
  param = unique_param,
  type_signal = unique_signal,
  smoothing = unique_smooth
)

method_comb$method_id <- seq(1,nrow(method_comb))

data <- left_join(df_all, method_comb, by = c("loc_method", "param", "type_signal", "smoothing"))

unique(is.na(data))

## Kalman smooth function -------------------------------------------------------
# run for each track, separating into 2-hour bins
kalman_track_smooth <- do.call("rbind", lapply(unique(data$method_id), function(method) {
  
  m_df <- data %>% filter(method_id == method)
  
  trial_loop <- do.call("rbind", lapply(unique(m_df$unique_trial_id), function(trial){
    
    result <- tryCatch({
    # select tag
    trial <- m_df %>% filter(unique_trial_id == trial)
    
    ## for x coordinates
    # make equal time intervals
    ts_x <- xts(trial$x.est, order.by = trial$test_local_time,
                tzone = "America/Bogota")
    
    # Make the time series regularly spaced with NAs for missing times
    r_ts_x <- seq(start(ts_x), end(ts_x), by = "1 sec")
    m_ts_x <- xts(merge(ts_x, xts(, order.by = r_ts_x)), order.by = r_ts_x)
    
    # fit state space model, with non-informative initial parameters
    ssm_x <- SSModel(m_ts_x$ts_x ~ SSMtrend(1, Q = list(matrix(NA), matrix(NA))), H = matrix(NA))
    fit2_x <- fitSSM(ssm_x, inits = c(0, 0, 0, 0), method = "BFGS")
    ks_result_x <- KFS(fit2_x$model, filtering = "state", smoothing = "state")
    
    # only get results for obs time stamps
    non_na_indices_x <- which(!is.na(m_ts_x$ts_x))
    
    k_smooth_results_x <- ks_result_x$alphahat
    k_smooth_sel_x <- k_smooth_results_x[non_na_indices_x]
    
    ## for y coordinates
    # make equal time intervals
    ts_y <- xts(trial$y.est, order.by = trial$test_local_time,
                tzone = "America/Bogota")
    
    # Make the time series regularly spaced with NAs for missing times
    r_ts_y <- seq(start(ts_y), end(ts_y), by = "1 sec")
    m_ts_y <- xts(merge(ts_y, xts(, order.by = r_ts_y)), order.by = r_ts_y)
    
    # fit state space model, with non-informative initial parameters
    ssm_y <- SSModel(m_ts_y$ts_y ~ SSMtrend(1, Q = list(matrix(NA), matrix(NA))), H = matrix(NA))
    fit2_y <- fitSSM(ssm_y, inits = c(0, 0, 0, 0), method = "BFGS")
    ks_result_y <- KFS(fit2_y$model, filtering = "state", smoothing = "state")
    
    # only get results for obs time stamps
    non_na_indices_y <- which(!is.na(m_ts_y$ts_y))
    
    k_smooth_results_y <- ks_result_y$alphahat
    k_smooth_sel_y <- k_smooth_results_y[non_na_indices_y]
    
    ## make dataframe
    data.frame(unique_trial_id = trial$unique_trial_id,
               reloc_id = trial$reloc_id,
               test_local_time = trial$test_local_time, 
               method_id = trial$method_id,
               x_ks = as.vector(k_smooth_sel_x),
               y_ks = as.vector(k_smooth_sel_y)) # smoothed estimate of states
  }, error = function(e){
    return(NULL)
  })
    return(result)
    
  }))
  
}))

kalman_data <- left_join(kalman_track_smooth, data, by = c("unique_trial_id", "reloc_id", "test_local_time", "method_id"))

# new coordinate estimate names
kalman_df <- kalman_data %>% select(-c(x.est, y.est, error_m)) %>% 
  rename(x.est = x_ks, y.est = y_ks)

# calculate difference between known and estimated locations
kalman_df$error_m <- raster::pointDistance(kalman_df[,c("x.est", "y.est")], kalman_df[,c("TestUTMx", "TestUTMy")], lonlat = F, allpairs = F)



## Cubic spline smooth -------------------------------------------------------

cubic_spline_smooth <- do.call("rbind", lapply(unique(data$method_id), function(method) {
  
  m_df <- data %>% filter(method_id == method)
  
  trial_loop <- do.call("rbind", lapply(unique(m_df$unique_trial_id), function(trial){
  
  track <- m_df %>% filter(unique_trial_id == trial)

    if (length(track$reloc_id) >= 5) {
      
      sp_x <- smooth.spline(x = track$test_local_time, y = track$x.est)
      sp_y <- smooth.spline(x = track$test_local_time, y = track$y.est)
      
      data.frame(unique_trial_id = track$unique_trial_id,
                 reloc_id = track$reloc_id,
                 method_id = track$method_id,
                 x_cs = sp_x[["y"]], # smoothed value of x
                 y_cs = sp_y[["y"]]) # smoothed value of y
      
    } else {
      
      data.frame(unique_trial_id = track$unique_trial_id,
                 reloc_id = track$reloc_id,
                 method_id = track$method_id,
                 x_cs = track$x.est, # smoothed value of x
                 y_cs = track$y.est) # smoothed value of y
      
    }
  }))
}))

cubic_spline_data <- left_join(cubic_spline_smooth, data, by = c("unique_trial_id", "reloc_id", "method_id"))

# new coordinate estimate names
cubic_spline_df <- cubic_spline_data %>% select(-c(x.est, y.est, error_m)) %>% 
  rename(x.est = x_cs, y.est = y_cs)

# calculate difference between known and estimated locations
cubic_spline_df$error_m <- raster::pointDistance(cubic_spline_df[,c("x.est", "y.est")], cubic_spline_df[,c("TestUTMx", "TestUTMy")], lonlat = F, allpairs = F)
