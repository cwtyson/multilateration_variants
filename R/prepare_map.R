## Prepare
library(tidyverse)

## Get grid points
source("./R/get_grid_points_fn.R")
source("./R/prepare_fn.R")
grid_points <- get_grid_points_fn(grid_points_folder = '/Users/tyson/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/My Drive/Zebby_tracking_field_data/grid_points',
                                        crs = 3308)

## Read in processed
dets <- readRDS("./data/dets_calibration_rounds_1-2.RDS") 

## Generate relevant combinations of arguments to ry
cbs <- data.frame(
  
  window = c(5,10,30),
  lag = c(0,-5,-15),
  dist_filter = c(175,325,1000), 
  func = c("mean","max","median")) %>% 
  expand.grid() %>% 
  filter(window > abs(lag)) %>% 
  mutate(window = paste(window, "secs"),
         lag = paste(lag, "secs"),
         func = as.character(func))

## Combination of arguments to run
args = list(cbs$window,
            cbs$lag,
            cbs$dist_filter,
            cbs$func)


## Prepare data using combinations of options
purrr::pmap(args,
            function(a,b,c,d) prepare_fn(dets = dets,
                                         grid_points = grid_points,
                                         tz = "UTC",
                                         crs = 3308,
                                         window = a,
                                         lag = b,
                                         dist_filter = c,
                                         sumFun = d))
