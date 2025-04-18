prepare_fn <- function(dets,
                       grid_points_folder,
                       tz,
                       crs,
                       window,
                       lag,
                       dist_filter,
                       sumFun){
  
  cat(window, "-", lag, "-", dist_filter, "-", sumFun, "\n")
  
  ## Function: 
  ## TODO - removes outliers 
  ## calculates running mean/max/median
  ## within rounded interval (window of rolling average) takes mean/max/median
  ## retains nodes within X meters of strongest node
  
  ## Apply sliding window and calculate specified summary stat
  fdets_prep <- dets %>%
    # slice(1:100) %>%
    dplyr::arrange(date_time) %>% 
    dplyr::group_by(tag, place, type, state, point, grid_point) %>% 
    dplyr::mutate(dt_r = lubridate::round_date(date_time,
                                               unit = window),
                  rssi_gp = runner::runner(x = .,
                                           k = window,
                                           lag = lag,
                                           idx = "date_time",
                                           f = function(x) match.fun(sumFun)(x$rssi),
                                           na_pad = FALSE)) %>% 
    dplyr::ungroup()
  
  ## Summarize based on rounded time
  fdets_prep_sum <- fdets_prep %>% 
    dplyr::group_by(tag,
                    place,
                    type,
                    state,
                    point,
                    dt_r,
                    grid_point)  %>%
    
    ## Summarize within rounded interval using specified summary function
    dplyr::summarise(rssi_gp_r = round(match.fun(sumFun)(rssi_gp, na.rm = T), 0),
                     .groups = "keep") %>% 
    
    ## Change groups
    dplyr::group_by(tag,
                    place,
                    type,
                    state,
                    point,
                    dt_r) %>% 
    
    ## Number of grid points in the interval
    dplyr::mutate(n_gp = n()) %>% 
    
    ## Arrange by dt_r, rssi_gp and then shuffle any ties
    dplyr::arrange(dt_r, desc(rssi_gp_r), runif(nrow(.))) %>% 
    
    ## Remove any intervals without at least 3 grid points
    dplyr::filter(n_gp >= 3) %>% 
    
    ## Join coordinates for grid points
    dplyr::left_join(grid_points, 
                     by = dplyr::join_by(grid_point)) %>% 
    
    ## Within each interval, distance from the first (strongest node) to the remainder
    dplyr::mutate(dist = sqrt((gp_x - dplyr::first(gp_x))^2 + (gp_y - dplyr::first(gp_y))^2)) %>% 
    
    ## Remove nodes farther away than specified filter
    dplyr::filter(dist <= dist_filter) %>% 
    
    ## Calculate remaining nodes within interval and remove intervals without at least 3
    dplyr::group_by(tag,
                    place,
                    type,
                    state,
                    point,
                    dt_r) %>% 
    dplyr::mutate(n_gp = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(n_gp >= 3) %>% 
    dplyr::as_tibble()  %>% 
    dplyr::mutate(sum_fun = sumFun,
                  window = window,
                  lag = lag,
                  dist_filter = dist_filter)  %>% 
    
    ## Group by:
    dplyr::group_by(tag, place, type, state, point, dt_r) %>% 
    
    ## Create interval id
    dplyr::mutate(int_id = dplyr::cur_group_id()) %>% 
    dplyr::arrange(int_id) %>% 
    dplyr::select(int_id, dplyr::everything())
  
  
  ## Save 
  saveRDS(fdets_prep_sum, paste0("./outputs/prepared/",sumFun,"_",window,"_",lag,"_",dist_filter,".RDS"))
  
}

