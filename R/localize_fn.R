localize_fn <- function(prepared_file,
                        rssi_dist_model_file,
                        reps = 100){
  
  cat("Starting to localize detection data\n")
  
  ## Read in prepared data
  dets_p <- readRDS(paste0("./outputs/prepared/", prepared_file))
  start_time = Sys.time()
  # rssi_dist_model_file = "/Users/tyson/Documents/academia/institutions/WUR/research/eswatini/eswatini_tracking/outputs/nls_mod.RDS"
  
  ## Read in model
  rssi_dist_mdl <- readRDS(rssi_dist_model_file) 
  
  ## Model table
  mdl_tab <- broom::augment(rssi_dist_mdl,
                            newdata = data.frame(rssi = seq(max(dets_p$rssi_gp_r), min(dets_p$rssi_gp_r), by = -1)),
                            se_fit = TRUE) %>% 
    dplyr::select(rssi_gp_r = rssi,
                  mean = .fitted,
                  sd = .se.fit) %>% 
    dplyr::distinct(rssi_gp_r,.keep_all = T)
  
  ## Process
  dets_p <- dets_p %>% 
    
    ## Join modeled rssi values
    left_join(mdl_tab, by = join_by(rssi_gp_r))
  
  ## Empty df
  tag_loc_est <- data.frame()
  
  ## Empty lists
  tag_loc_est_list <- list()
  cov_list <- list()
  pe_df_list <- list()
  
  cat("\n Intervals to localize:", max(dets_p$int_id),  "- with -", reps, "resamples each") 
  
  ## Set progress bar
  pb <- txtProgressBar(min = 0, max = max(dets_p$int_id), style = 3)
  
  ## For each interval, localize X times
  for(int in unique(dets_p$int_id)){
    
    ## Progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, which(int == 1:unique(max(dets_p$int_id))))
    
    ## Subset detections
    dets_p_int = dets_p[dets_p$int_id == int,]
    
    ## Point estimates
    pe <- list()
    
    for(i in 1:reps){
      
      ## Sample within interval
      dets_p_int <- dets_p_int %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(dist_est_samp = round(10^(rnorm(n = 1, 
                                                      mean = mean,
                                                      sd = sd)),
                                            0)) %>%
        dplyr::ungroup()
      
      tryCatch(
        expr = {
          
          # Non-linear test to optimize the location estimate given pairwise distances
          nls_mod <- suppressWarnings(nls(dist_est_samp ~ raster::pointDistance(data.frame(gp_x, gp_y), 
                                                                                c(x_solution, y_solution),
                                                                                lonlat = F, 
                                                                                allpairs = T), 
                                          data = dets_p_int,
                                          start = list(x_solution = dets_p_int[1,]$gp_x,
                                                       y_solution = dets_p_int[1,]$gp_y),
                                          control = nls.control(warnOnly = T,
                                                                maxiter=20,
                                                                minFactor=1/3000)))
          
          ## Combine point estimates
          pe[[i]] <- c(coef(nls_mod)[1],
                       coef(nls_mod)[2])
          
        },
        
        ## Print error message
        error = function(e)  {cat("ERROR :",conditionMessage(e), "\n")})
      
    }
    
    
    ## Combine point estimates
    pe_df <- do.call(rbind,pe) %>% 
      data.frame()
    
    ## Ellipse info from X repeated multilateration estimates
    ell.info <- cov.wt(cbind(pe_df$x_solution, pe_df$y_solution))
    
    ## Center
    x_est <- ell.info$center[1]
    y_est <- ell.info$center[2]

    ## Axes
    eigen.info <- eigen(ell.info$cov)
    lengths <- sqrt(eigen.info$values * 2 * qf(.95, 2, reps-1))
    a = lengths[1]
    b = lengths[2]
    
    ## Orientation
    u <- eigen.info$vectors[,1] # major axis eigenvector
    theta <- atan2(u[2],u[1]) # angle from x-axis in radians
    theta <- theta * 360/(2*pi) # angle from x-axis in degrees
    theta <- theta-90 # angle from y-axis in degrees
    
    ## Combine estimated locations with summary information
    tag_int_loc_est <- dets_p_int %>% 
      select(where(~ n_distinct(.) == n_distinct(dets_p_int$int_id))) %>%     
      distinct() %>% 
      mutate(mean_RSSI = mean(dets_p_int$rssi_gp_r),
             max_RSSI = max(dets_p_int$rssi_gp_r),
             x_est = x_est,
             y_est = y_est,
             semi_major = a,
             semi_minor = b,
             orientation = theta)
    
    ##  Get raw points
    pe_df_p <- tag_int_loc_est %>% 
      select(int_id,
             tag,
             place,
             type,
             point,
             state)  %>% 
      left_join(pe_df %>% 
                  mutate(int_id = tag_int_loc_est$int_id),
                by = join_by(int_id))
    
    ## Add elements to list
    tag_loc_est_list[int] <- list(tag_int_loc_est)
    pe_df_list[int] <- list(pe_df_p)
    cov_list[int] <- list(ell.info$cov)
    
  }
  
  ## Combine estimated locations, raw point estimates, and covariance matrices
  tag_loc_est <- list(tag_loc_est_list %>% 
                        do.call(rbind,.),
                      pe_df_list %>% 
                        do.call(rbind,.),
                      cov_list %>% 
                        do.call(rbind,.))
  
  # ggplot() +
  #   geom_point(aes(x=x_solution,
  #                  y=y_solution),,
  #              data= tag_loc_est[[2]]) +
  #   ggforce::geom_ellipse(aes(x0=x_est,
  #                             y0=y_est,
  #                             a=semi_major,
  #                             b=semi_minor,
  #                             angle=orientation),
  #                         data = tag_loc_est[[1]],
  #                         color = "red") +
  #   geom_point(aes(gp_x,gp_y),data=dets_p %>%
  #                select(gp_x,gp_y) %>%
  #                distinct()) +
  #   facet_wrap(int_id~.) +
  #   theme_minimal()
  
  cat("\n Finished localizing after:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "minutes")
  
  ## End progress bar for all tags
  close(pb)
  
  ## Save RDS
  saveRDS(tag_loc_est, 
          paste0("./outputs/localized/",prepared_file))
}
