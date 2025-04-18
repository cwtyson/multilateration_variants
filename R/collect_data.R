## Get detections from validation walks

## Housekeeping
library(tidyverse)
library(janitor)
library(sf)
library(fuzzyjoin)
library(stringr)

## Get raw detections ######

## Tag logs to define dates for filtering
log_files <- list.files("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/",
                        pattern = "point_logs",
                        full.names = T,
                        recursive = TRUE) 

log_files <- stringr::str_subset(log_files, fixed('~$'), negate = TRUE)

## Reformat
logs <- do.call(bind_rows, lapply(log_files, function(x) readxl::read_excel(x,col_types="text")) ) %>% 
  mutate(start_time = dmy_hm(paste(date, start_time),tz="Australia/Broken_Hill"),
         end_time = dmy_hm(paste(date, end_time),tz="Australia/Broken_Hill"))

## Times - convert to UTC and adjust by one hour on either side
min_time <- with_tz(min(logs$start_time,na.rm=T), "UTC") - 3600
max_time <- with_tz(max(logs$end_time,na.rm=T),"UTC") + 3600

## Connect to data base back end
conn <- DBI::dbConnect(duckdb::duckdb((dir = "/Users/tyson/Documents/academia/institutions/WUR/research/CTT_data/Zebra Finches in Australia/zebby_database.duckdb")),
                       read_only= T)

## Tags from both rounds
tags1 <- readxl::read_excel("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round1_07_Sept_2024/walking_validation_07Sept2024_tag_log.xlsx") %>% 
  pull(tag) 

tags2 <- readxl::read_excel("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round2_13_Sept_2024/walking_validation_13_14_Sept2024_tag_log.xlsx") %>% 
  pull(tag) 

tags <- c(tags1,tags2)

## Read in detections from postgres database
dets_raw <- dplyr::tbl(conn, "raw") %>%
  
  ## Keep only station ids matching the specified filter
  dplyr::filter(station_id %in% c("D82AA0A12259", "4BA80216EAEB")) %>%
  
  ## Keep tags
  dplyr::filter(tag_id %in% tags) %>%
  
  ## Filter based on time
  dplyr::filter(time > min_time) %>%
  dplyr::filter(time < max_time) %>%
  
  ## Collect filitered data from database
  dplyr::collect() %>%
  
  ## Distinct values
  dplyr::distinct(tag_id,
                  node_id,
                  time,
                  .keep_all = T) %>%
  ## Select and rename
  dplyr::transmute(node = toupper(node_id),
                   date_time = lubridate::with_tz(time, tz = "Australia/Broken_Hill"),
                   tag = tag_id,
                   rssi = tag_rssi) %>%
  arrange(date_time)

## Filter
dets_f <- dets_raw %>%
  
  ## Remove impossible RSSI values
  filter(rssi < 0)

## Node codes
node_codes <- readxl::read_excel("/Users/tyson/Google Drive/My Drive/Zebby_tracking_field_data/nodes/node_codes_20230906.xlsx") %>% 
  clean_names() %>% 
  transmute(node = toupper(node),
            node_number = as.character(node_number))

## Node deployment log
node_log <- readxl::read_excel("/Users/tyson/Google Drive/My Drive/Zebby_tracking_field_data/nodes/node_deployment_log_20240927.xlsx") %>% 
  dplyr::transmute(grid_point,
                   node_number = as.character(round(as.numeric(node_number),0)),
                   deployment_time = lubridate::dmy_hm(paste(start_date, start_time), 
                                                       tz = "Australia/Broken_Hill"),
                   removal_time = lubridate::dmy_hm(paste(end_date, end_time),
                                                    tz = "Australia/Broken_Hill")) %>% 
  ## Join node node
  dplyr::left_join(node_codes,
                   by  = "node_number") %>% 
  
  dplyr::select(node,
                grid_point,
                node_number,
                date_time = deployment_time,
                removal_time)  

## Convert to data.table
nodes <- data.table::data.table(node_log, key = c("node", "date_time"))
dets_f <- data.table::data.table(dets_f, key = c("node", "date_time"))

## Rolling join node log to node records
dets_f <- nodes[dets_f, roll = Inf]

## Remove detections without a grid point and where date time is after removal time
dets_f1 <- dets_f %>%
  dplyr::filter(!is.na(grid_point),
                (date_time < removal_time | is.na(removal_time))) %>%
  
  ## Arrange and select
  dplyr::arrange(tag,
                 date_time) %>%
  dplyr::select(grid_point,
                tag,
                date_time,
                rssi) %>%
  data.frame() 

saveRDS(dets_f1,
        "./data/dets_raw.RDS")

## Collect from round 1 ########

## Read in data
dets_f1 <- readRDS("./data/dets_raw.RDS")

## Tags
tag_log <- readxl::read_excel("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round1_07_Sept_2024/walking_validation_07Sept2024_tag_log.xlsx")

## Just tags
tags <- tag_log %>% 
  pull(tag)  

## Keep only legitimate tags
dets_r2  <- dets_f1 %>% 
  filter(tag %in% tags)

## Point log
field_log_intervals <- readxl::read_xlsx("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round1_07_Sept_2024/stopping_point_logs.xlsx") %>% 
  mutate(start_time = dmy_hm(paste(date, start_time),tz="Australia/Broken_Hill")+60,
         end_time = dmy_hm(paste(date, end_time),tz="Australia/Broken_Hill"),
         interval = lubridate::interval(start_time, end_time)) %>% 
  select(name,point,start_time,end_time) %>% 
  mutate(state = "S")

int_log <- interval(min(field_log_intervals$start_time),
                    max(field_log_intervals$end_time))

## Add tag info
dets_j <- dets_r2 %>% 
  left_join(tag_log)

## Empty data frame
field_cal_dets_j_r1 <- data.frame()
for(p in unique(dets_j$name)){
  
  cat(p,"\n")
  
  ## Subset name  
  field_cal_dets_p <- dets_j[dets_j$name == p,]
  
  field_log_interval_p <- field_log_intervals[field_log_intervals$name == p,]
  
  ## Fuzzy join based on dates
  dets_fuz_p <- field_cal_dets_p %>% 
    fuzzy_left_join(field_log_interval_p,
                    by = c(
                      "name" = "name",
                      "date_time" = "start_time",
                      "date_time" = "end_time"
                    ),
                    match_fun = list(`==`, `>=`, `<=`)) %>% 
    select(tag,date_time,gp=grid_point,rssi,name = name.x,place,type,point,state) %>% 
    mutate(state = ifelse(is.na(state),"M",state))
  
  ## Combine
  field_cal_dets_j_r1 <- dplyr::bind_rows(field_cal_dets_j_r1, dets_fuz_p)
  
}

## Keep dets between start and end
fcdf_r1 <- field_cal_dets_j_r1 %>% 
  filter(date_time > min(field_log_intervals$start_time)) %>% 
  filter(date_time < max(field_log_intervals$end_time)) 


## Collect from round 2 ########

## Tags
tag_log <- readxl::read_excel("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round2_13_Sept_2024/walking_validation_13_14_Sept2024_tag_log.xlsx")

## Just tags
tags <- tag_log %>% 
  pull(tag)  

## Keep only legitimate tags
dets_r2  <- dets_f1 %>% 
  filter(tag %in% tags)

## Point log -- Check how Noelle and Tis did this
field_log_intervals <- readxl::read_xlsx("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round2_13_Sept_2024/random_stopping_point_logs.xlsx") %>% 
  mutate(start_time = dmy_hm(paste(date, start_time),tz="Australia/Broken_Hill"),
         end_time = dmy_hm(paste(date, end_time),tz="Australia/Broken_Hill"),
         interval = lubridate::interval(start_time, end_time)) %>% 
  select(name,point,start_time,end_time,state=type)

int_log <- interval(min(field_log_intervals$start_time),
                    max(field_log_intervals$end_time))

## Add tag info
dets_j <- dets_r2 %>% 
  left_join(tag_log)

## Empty data frame
field_cal_dets_j_r2 <- data.frame()
for(p in unique(dets_j$name)){
  
  cat(p,"\n")
  
  ## Subset name  
  field_cal_dets_p <- dets_j[dets_j$name == p,]
  
  field_log_interval_p <- field_log_intervals[field_log_intervals$name == p,]
  
  ## Fuzzy join based on dates
  dets_fuz_p <- field_cal_dets_p %>% 
    fuzzy_left_join(field_log_interval_p,
                    by = c(
                      "name" = "name",
                      "date_time" = "start_time",
                      "date_time" = "end_time"
                    ),
                    match_fun = list(`==`, `>=`, `<=`)) %>% 
    select(tag,date_time,gp=grid_point,rssi,name = name.x,place,type,point,state) %>% 
    mutate(state = ifelse(is.na(state),"M",state))
  
  ## Combine
  field_cal_dets_j_r2 <- dplyr::bind_rows(field_cal_dets_j_r2, dets_fuz_p)
  
}

## Keep dets between start and end
fcdf_r2 <- field_cal_dets_j_r2 %>% 
  filter(date_time > min(field_log_intervals$start_time)) %>% 
  filter(date_time < max(field_log_intervals$end_time))

## Combine rounds
dets_all <- bind_rows(fcdf_r1 %>% 
                        mutate(round = "one"),
                      fcdf_r2 %>% 
                        mutate(round = "two")) %>% 
  
  ## Reorder columns
  select(point,name,tag,place,type,state, date_time,grid_point = gp,rssi)
  
## Save
saveRDS(dets_all, 
        "./data/dets_calibration_rounds_1-2.RDS")

## Check detections ######

## Detections
dets_all <- readRDS("./data/dets_calibration_rounds_1-2.RDS")

## Processed detecionts
dets_all_p <- dets_all %>% 
  group_by(round, point) %>% 
  mutate(seconds = max(date_time) - min(date_time)) %>% 
  na.omit() %>% 
  filter(type == "hybrid") %>% 
  filter(round == "two") %>% 
  group_by(point,type,place,tag) %>% 
  mutate(dets = n(),
         dets_rate = n()/as.numeric(seconds),
         mean_rssi = mean(rssi),
         max_rssi = max(rssi),
         nodes = n_distinct(gp)) %>% 
  distinct(point,type,place,tag,.keep_all = T) %>% 
  select(point,tag,type,place,seconds,dets:nodes)

dets_all_pl <- dets_all_p %>% 
  pivot_longer(cols = c(dets_rate:nodes))

## RSSI
ggplot(dets_all_pl) +
  geom_jitter(aes(x=place,
                  y=value,
                  color=tag)) +
  facet_wrap(.~name,
             scales="free") +
  theme_minimal()


ggplot(dets_all_p) +
  geom_jitter(aes(x=tag,y=dets,color=max_rssi)) +
  facet_grid(.~place,scales="free") +
  theme_minimal()



## Paths 
wv_routes1 <- sf::st_read("./data/2024/walking_validation/round1_07_Sept_2024/wv_field_tracks.GPX",
                          layer = "track_points") %>% 
  dplyr::transmute(name = dplyr::case_when(track_fid == 0 ~ "noelle",
                                           track_fid == 1 ~ "tis",
                                           track_fid == 2 ~ "chris"),
                   date_time = with_tz(time, tz="Australia/Broken_Hill")) %>% 
  group_by(name) %>% 
  mutate(speed = as.numeric(st_distance(geometry, lead(geometry), by_element = T)/5),
         point = 1:n())

## Paths 
wv_routes2 <- sf::st_read("./data/2024/walking_validation/round2_13_Sept_2024/wv_field_tracks.GPX",
                          layer = "track_points") %>% 
  dplyr::transmute(name = dplyr::case_when(track_fid == 0 ~ "noelle",
                                           track_fid == 1 ~ "tis",
                                           track_fid == 2 ~ "chris"),
                   date_time = with_tz(time, tz="Australia/Broken_Hill")) %>% 
  group_by(name) %>% 
  mutate(speed = as.numeric(st_distance(geometry, lead(geometry), by_element = T)/5),
         point = 1:n()) 

## Walked paths
walked <- bind_rows(wv_routes1 %>% 
                      mutate(round = "one"),
                    wv_routes2 %>% 
                      mutate(round = "two"))


## Point log
field_log_intervals1 <- readxl::read_xlsx("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round1_07_Sept_2024/stopping_point_logs.xlsx") %>% 
  mutate(start_time = dmy_hm(paste(date, start_time),tz="Australia/Broken_Hill"),
         end_time = dmy_hm(paste(date, end_time),tz="Australia/Broken_Hill"),
         interval = lubridate::interval(start_time, end_time)) %>% 
  mutate(state = "S",
         round="one") %>% 
  select(round,name,route,point,state,start_time,end_time) %>% 
  
  ## Tis' watch was ~25s fast, so cutoff 30s from her time
  mutate(end_time = case_when(name == "Tis" ~ end_time -30,
                              TRUE ~ end_time))


## Point log -- Check how Noelle and Tis did this
field_log_intervals2 <- readxl::read_xlsx("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/round2_13_Sept_2024/random_stopping_point_logs.xlsx") %>% 
  mutate(start_time = dmy_hm(paste(date, start_time),tz="Australia/Broken_Hill"),
         end_time = dmy_hm(paste(date, end_time),tz="Australia/Broken_Hill"),
         interval = lubridate::interval(start_time, end_time),
         route = str_split_i(point,"_",1),
         round="two") %>% 
  select(round,name,route,point,start_time,end_time,state=type)

## Field log routes
route_log <- bind_rows(field_log_intervals1,field_log_intervals2) %>% 
  group_by(round, name,route) %>% 
  summarise(start_time = min(start_time),
            end_time = max(end_time)) %>% 
  arrange(start_time) %>% 
  mutate(name = tolower(name))


## Add detection information - rounding GPS points and detections to 5 second interval and joining
walked_p <- walked %>% 
  mutate(dt_r = round_date(date_time, unit = "5 secs")) %>% 
  distinct(name,dt_r,.keep_all = T) %>% 
  arrange(dt_r) %>%
  group_by(dt_r) %>% 
  mutate(pt_time = cur_group_id()) %>% 
  mutate(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2])

## Join route times
walked_p2 <- walked_p %>%
  # st_drop_geometry() %>%
  ungroup() %>% 
  select(round,name,date_time, x,y) %>% 
  fuzzy_left_join(route_log,
                  by = c(
                    "round" = "round",
                    "name" = "name",
                    "date_time" = "start_time",
                    "date_time" = "end_time"
                  ),
                  match_fun = list(`==`,`==`, `>=`, `<=`)) %>% 
  na.omit() %>% 
  select(round = round.x,
         route,
         x,
         y)

## Read in processed
dets <- readRDS("./data/2024/walking_validation/dets_rounds1-2.RDS")

## Plot
ggplot(dets) +
  geom_point(aes(x=date_time,
                 y=tag)) +
  facet_grid(place+state+type~.,scales="free")

## At each calibration point, get strongest node
dets_cp <- dets %>% 
  na.omit() %>% 
  
  # ## Remove weak detections
  # filter(rssi > -90) %>%
  # 
  group_by(round,place,type,point,state,gp,tag) %>% 
  summarise(rssi = round(max(rssi)),
            dets = n()) %>% 
  group_by(round,point,place,type,state,gp,tag) %>% 
  arrange(round,point,place,type,state,desc(rssi)) %>% 
  group_by(round,point,state,tag) %>% 
  slice_head() %>% 
  data.frame()  %>% 
  
  ## Tags per point
  group_by(round,point,state,gp) %>% 
  mutate(tags = n())

## Points
gps_pt_files <- list.files("./data/2024/walking_validation",recursive = T,pattern="points",full.names = T)

## Calibration points
gps_pts <- do.call(bind_rows, lapply(gps_pt_files, function(x) st_read(x) %>% 
                                       janitor::clean_names()) ) %>% 
  select(point=name)

## Grid points
gps <- st_read("/Users/tyson/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/My Drive/Zebby_tracking_field_data/grid_points/grid_points_20240916.GPX") %>% 
  select(grid_point = name)  %>% 
  filter(grepl("Gp",grid_point)) %>% 
  mutate(grid_point = gsub("Gp ", "", grid_point))  %>% 
  mutate(gp_x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         gp_y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2])

## Join coordinates of grid point to calibration points strongest detecting node
dets_all <- dets %>% 
  left_join(gps_pts) %>% 
  mutate(cp_x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         cp_y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
  left_join(gps, by = c("gp" = "grid_point")) %>% 
  select(round,
         point,
         tag,
         cp_x,
         cp_y,
         gp,
         gp_x,
         gp_y
  ) %>% 
  na.omit()

## Visualize detections at each point
ggplot(dets_all ) +
  
  # geom_point(aes(x=x,y=y,color=route),
  #            data = walked_p2) +
  # 
  
  geom_point(aes(x=gp_x,
                 y=gp_y),
             color="blue",
             alpha=0.2,
             data = gps) +
  
  geom_segment(aes(x=cp_x,
                   y=cp_y,
                   xend=gp_x,
                   yend=gp_y)) +
  # geom_text(aes(x=cp_x,
  #               y=cp_y,
  #               label=point),
  #           color = "red") +
  # 
  
  facet_wrap(gp~.) +
  theme_void() 


## Join coordinates of grid point to calibration points strongest detecting node
dets_cp_p <- dets_cp %>% 
  left_join(gps_pts) %>% 
  mutate(cp_x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         cp_y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
  left_join(gps, by = c("gp" = "grid_point")) %>% 
  
  ## Distance between points
  mutate(
    dist = geosphere::distHaversine(cbind(gp_x, gp_y), cbind(cp_x, cp_y))
  )


ggplot(dets_cp_p ) +
  
  geom_text(aes(x=gp_x,
                y=gp_y,
                label=grid_point),
            data = gps) +
  
  # geom_point(aes(x=gp_x,
  #                y=gp_y),
  #            size = ) +
  
  geom_point(aes(x=x,y=y,color=route),
             data = walked_p2) +
  
  geom_segment(aes(x=cp_x,
                   y=cp_y,
                   xend=gp_x,
                   yend=gp_y)) +
  geom_text(aes(x=cp_x,
                y=cp_y,
                label=point),
            color = "red") +
  theme_minimal() +
  facet_grid(~round)

## Distance / rssi
ggplot(dets_cp_p %>% 
         filter(round == "one")) +
  geom_text(aes(x=dist,
                y=max,
                label=point)) +
  facet_grid(type~place)

## Plot features of RSSI in different states ######

## State key:
# S - stationary, pole was held at the point motionless
# A - active, pole was held at the point and rotated continuously
# M - moving, moving between points

## Read in data
dets_all_p <- readRDS("./data/2024/walking_validation/dets_rounds1-2.RDS") %>% 
  mutate(dt_r = lubridate::floor_date(date_time, unit = "1min"),
         point = ifelse(is.na(point),"not",point)) %>% 
  group_by(tag, state, dt_r) %>% 
  mutate(group_id = cur_group_id(),
         count = n()) %>% 
  select(round,tag,place,type,state,point,group_id,dt_r,gp,rssi)

## Number of topp nodes to keep
top_node_num = 5

