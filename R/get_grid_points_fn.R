get_grid_points_fn <- function(grid_points_folder,
                                     crs){
  
  
  ## Get all grid point files
  all_grid_points <- sort(list.files(grid_points_folder,
                                     full.names = TRUE,
                                     pattern = "points"),
                          decreasing = TRUE)
  
  ## Read in each file and combine
  grid_points <- data.frame()
  for(gp_f in all_grid_points){
    
    ## Get grid point coordinates
    grid_points_f <- suppressWarnings(sf::read_sf(gp_f) %>% 
                                        sf::st_transform(crs) %>% 
                                        janitor::clean_names() %>% 
                                        dplyr::transmute(grid_point = name,
                                                         gp_x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
                                                         gp_y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
                                        sf::st_drop_geometry())
    
    grid_points <- bind_rows(grid_points,
                             grid_points_f) 
    
  }
  
  ## Keep unique points
  grid_points <- grid_points %>% 
    distinct(grid_point,
             .keep_all = T)
  
  ## Rename grid points by removing "Gp "
  grid_points$grid_point <- gsub(pattern = "Gp ",
                                 replacement = "",
                                 grid_points$grid_point)
  
  
  return(grid_points)
}

