## Transform downloaded RSS data to x and y coordinates

library(tidyverse)
library(lubridate)
library(sf)
library(cowplot)

# define select function to avoid confusion
select <- dplyr::select

## 0. Data prep ----------------------------------------------------------------

# import
beep_files <- list.files("beep_data/")

beep_list <- lapply(beep_files, function(file) {
  read.csv(paste0("beep_data/", file))
})

beep_df <- bind_rows(beep_list)

# remove RadioId to avoid duplicates
df <- beep_df %>% select(-RadioId) %>%
  distinct()

# check NAs in data, remove
unique(is.na(df))
df <- df %>% filter(!is.na(Time.local))

# convert time
# first format into POSIXct
df$Time.format <- as.POSIXct(strptime(df$Time.local,
                                      format = "%Y-%m-%dT%H:%M:%SZ", 
                                      tz = "GMT"))

# local time zone
df$Time.local <- with_tz(df$Time.format, tzone = "America/Bogota")

# check tz
df_check <- df %>% mutate(h = hour(Time.local))

(tz_check <- ggplot() +
    geom_density(data = df_check, aes(x = h)))

## add relevant tag info
tags <- read.csv("deployed_tags.csv")

df <- left_join(df, tags, by = "TagId")

# remove any detections from before attachment
df$Time.tag.attachment_format <- paste0(df$attachment_date, " ", df$attachment_hour, ":00")
df$Time.tag.attachment <- as.POSIXct(strptime(df$Time.tag.attachment_format,
                                              format = "%m/%d/%Y %H:%M:%S",
                                              tz = "America/Bogota"))

df <- df %>% filter(Time.local > Time.tag.attachment) %>% select(-Time.tag.attachment_format)

## 1 . Signal selection --------------------------------------------------------
# Time window of interest must be set by researcher,
# and also depends on frequency of signal emission

# In our case, Lifetags are set to emit every 2 s, Powertags every 60 s

# First we will explore how time difference between signals is distributed 
# in the data

# calculate time difference in node detections
dt_times <- df %>% 
  group_by(TagId, NodeId) %>% # for each ind and node
  arrange(Time.local) %>% # order by time
  mutate(Time.dif = Time.local - lag(Time.local)) %>% # calculate dt
  select(TagId, NodeId, Time.local, tag_type, Time.dif)

# NA times = 0
dt_times[is.na(dt_times$Time.dif),]$Time.dif <- as.duration(seconds(0))
unique(is.na(dt_times))

# Vis dt
quantile(dt_times$Time.dif, probs = 0.95)

(dt_g <- ggplot() +
    geom_density(data = dt_times, aes(x = Time.dif, color = tag_type)) +
    xlim(0,137) + # zoom into 95% values
    facet_wrap(~TagId, scales = "free") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()))

# distribution by tag type
power_dt_times <- dt_times %>% filter(tag_type == "PowerTag")
life_dt_times <- dt_times %>% filter(tag_type == "LifeTag")

quantile(power_dt_times$Time.dif, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.95, 1))
quantile(life_dt_times$Time.dif, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.95, 1))


## Now we aggregate time into intervals

# identify start times for each tag
tag_start_times <- df %>% 
  group_by(TagId) %>% # for each ind and node
  slice(which.min(Time.local)) %>%
  mutate(Tag.start.time = Time.local) %>%
  select(TagId, Tag.start.time)

# add aggregated time intervals
df_time <- left_join(df, tag_start_times, by = c("TagId"))


df_time <- df_time %>% mutate(Tag.start.time.diff = as.numeric(Time.local - Tag.start.time),
                              reloc.time.group = case_when(tag_type == "PowerTag" ~ floor(Tag.start.time.diff/120), # 120 s
                                                        tag_type == "LifeTag" ~ floor(Tag.start.time.diff/10)), # 10 s
                              bin.2h = floor(Tag.start.time.diff/(60*60*2))) # Create column for 2h bins for smoothing below

# choose maximum by node within time interval
max_df <- df_time %>% group_by(TagId, NodeId, reloc.time.group) %>%
  summarise(Time.local = min(Time.local), # time as the start of relocation interval
            maxRSS = max(TagRSSI),
            avgRSS = mean(TagRSSI),
            sdRSS = sd(TagRSSI),
            n.det = n()) %>%
  ungroup()

quantile(max_df$n.det)

# create unique relocation ids
reloc_df <- df_time %>% group_by(TagId, reloc.time.group) %>%
  summarise(reloc = cur_group_id()) %>%
  ungroup()

# join to calculated max
max_df <- left_join(max_df, reloc_df, by = c("TagId", "reloc.time.group"))

# complement with data
j_df <- df_time %>% select(TagId, NodeId, tag_type, tag_attachment_method, 
                           attachment_date, species, id, attachment_hour, 
                           Tag.start.time, reloc.time.group, bin.2h) %>% distinct()

sel_df <- left_join(max_df, j_df, by = c("TagId", "NodeId", "reloc.time.group"))

unique(is.na(sel_df))

## 2. Smoothing ----------------------------------------------------------------

# separate unique tags
unique_tags <- unique(sel_df$TagId)

# run by node for each track, separating into 2-hour bins
cubic_smoother <- do.call("rbind", lapply(unique_tags, function(tag) {
  
  # select tag
  track <- sel_df %>% filter(TagId == tag)
  
  bin_loop <- do.call("rbind", lapply(unique(track$bin.2h), function(bin){
    
    # select 2h bin
    bin_df <- track %>% filter(bin.2h == bin)
    
    node_loop <- do.call("rbind", lapply(unique(bin_df$NodeId), function(node) {
      
      # select node
      node_seq <- bin_df %>% filter(NodeId == node)
      
      if (nrow(node_seq) >= 5) {
        
        sp <- smooth.spline(x = node_seq$Time.local, y = node_seq$maxRSS, cv = TRUE)
        
        data.frame(TagId = node_seq$TagId,
                   NodeId = node_seq$NodeId,
                   Time.local = node_seq$Time.local,
                   maxRSS_cs = sp[["y"]]) # smoothed value
        
      } else {
        # if less than 5 consecutive relocs, have to keep original read
        data.frame(TagId = node_seq$TagId,
                   NodeId = node_seq$NodeId,
                   Time.local = node_seq$Time.local,
                   maxRSS_cs = node_seq$maxRSS) # smoothed value
        
      }
      
    }))
  }))
}))


# bind result to dataframe
cubic_results <- left_join(sel_df, cubic_smoother, by = c("TagId", "NodeId", "Time.local"))

# check join worked ok
unique(is.na(cubic_results))

# save
write_csv(cubic_results, "smoothed_tracks/cs_tracks.csv")

## 3. RSS to distance ----------------------------------------------------------------

# convert maxRSS to distance with 
# general exponential decay model calculated with maxRSS and Cubic smoothed values before

a <- 28.62524
S <- 0.008464565
K <- -102.4698

gen_df <- cubic_results %>% mutate(e_dist = (log(maxRSS_cs - K) - log(a)) / -S)

# remove NAs from df, when signal was lower than K
gen_filt_df <- gen_df %>% filter(!is.na(e_dist))

# lost one tag here, check which one
base::setdiff(gen_df$TagId, gen_filt_df$TagId)

no_signals_tag <- cubic_results %>% filter(TagId == "4C4B1E2D")
range(no_signals_tag$maxRSS)
range(no_signals_tag$maxRSS_cs)

# negative distances are assumed to be 0 (too close to node)
gen_relocs <- gen_filt_df %>% mutate(e_dist = case_when(e_dist < 0 ~ 0,
                                                     e_dist >=0 ~ e_dist))

# cannot use relocations with less than 3 nodes
n_nodes_relocs <- gen_relocs %>% group_by(reloc) %>% 
  summarise(n_nodes = n()) %>%
  ungroup()

relocs_over3_nodes <- n_nodes_relocs[n_nodes_relocs$n_nodes >2,]$reloc

relocs_under3_nodes <- n_nodes_relocs[n_nodes_relocs$n_nodes <= 3,]$reloc # keep

gen_relocs <- gen_relocs %>% filter(reloc %in% relocs_over3_nodes)

## total 
# n signals = 974791 (length(gen_df$Time.local)) n relocs = 120989 (length(unique(gen_df$reloc)))
## above h asymptote (K)
# n signals = 316974 length(gen_filt_df$Time.local) n relocs = 100483 length(unique(gen_filt_df$reloc))
## at least 3 nodes detected
#  n relocs = 43774 length(relocs_over3_nodes)

## compare with model varying parameters by node -------------------------------

# import parameter values by node
node_param <- read.csv("max_rss_cs_by_node.csv")
node_param <- node_param %>% select("NodeId", "a", "S", "K") %>% distinct()

node_df <- cubic_results

node_df <- left_join(node_df, node_param, by = "NodeId")
node_df[is.na(node_df$a),]$a <- 28.62524
node_df[is.na(node_df$S),]$S <- 0.008464565
node_df[is.na(node_df$K),]$K <- -102.4698

unique(is.na(node_df))

node_df$e_dist <- (log(node_df$maxRSS_cs - node_df$K) - log(node_df$a)) / -node_df$S

# remove NAs from df, when signal was lower than K
node_df_filt <- node_df %>% filter(!is.na(e_dist))

# lost one tag here, check which one
base::setdiff(node_df$TagId, node_df_filt$TagId)

# negative distances are assumed to be 0 (too close to node)
node_relocs <- node_df_filt %>% mutate(e_dist = case_when(e_dist < 0 ~ 0,
                                                        e_dist >=0 ~ e_dist))

# cannot use relocations with less than 3 nodes
n_nodes_relocs_n <- node_relocs %>% group_by(reloc) %>% 
  summarise(n_nodes = n()) %>%
  ungroup()

relocs_over3_nodes_n <- n_nodes_relocs_n[n_nodes_relocs_n$n_nodes >2,]$reloc

relocs_under3_nodes_n <- n_nodes_relocs_n[n_nodes_relocs_n$n_nodes <= 3,]$reloc # keep

node_relocs <- node_relocs %>% filter(reloc %in% relocs_over3_nodes_n)

## total 
# n signals = 974791 (length(node_df$Time.local)) n relocs = 120989 (length(unique(node_df$reloc)))
## above h asymptote (K)
# n signals = 390915 length(node_df_filt$Time.local) n relocs = 103954length(unique(node_df_filt$reloc))
## at least 3 nodes detected
#  n relocs = 56840 length(relocs_over3_nodes_n)

## vis
(rss_d <- ggplot() +
    geom_vline(data = node_param, aes(xintercept = K), color = "#454545", linetype = "dashed") +
    geom_vline(xintercept = -102.4698, lwd = 1, color = "blue") + # horizontal asymptote
    geom_density(data = cubic_results, aes(x = maxRSS), color = "#454545") +
    geom_density(data = cubic_results, aes(x = maxRSS_cs), color = "#FF1493") +
    facet_wrap(~TagId, scales = "free", ncol = 2) +
    labs(x = "maxRSS (dB)", y = "Density") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()))


## 4. Localization ----------------------------------------------------------------

## add node UTMS
Nodes <- read.csv("Nodes.csv")

Nodes <- Nodes %>% select(-c(Own_NodeId, Latitude, Longitude)) # remove unnecessary columns

# careful, some node ids revert to different names when ends in e followed by number
Nodes[Nodes$NodeId == "3288E6",]$NodeId <- "3288e6"

# join to relocs df
node_relocs <- left_join(node_relocs, Nodes, by = "NodeId")

unique(is.na(node_relocs))

# vector of unique relocs
unique_relocations <- unique(node_relocs$reloc)

# run trilateration with 3 strongest rss nodes
trilat_rssi <- do.call("rbind", lapply(unique_relocations, function(reloc_id) {
  
  # Separate by relocations
  reloc_single <- node_relocs %>% filter(reloc == reloc_id)
  
  # Get only 3 top nodes
  reloc_top_rssi <- reloc_single %>% slice_max(maxRSS_cs, n = 3, with_ties = FALSE) # only 3, no ties
  
  # Get location of nodes with signal for relocation
  node_UTM <- reloc_top_rssi %>% select(NodeUTMx, NodeUTMy) %>% distinct()
  
  # Calculate no nodes for the test
  no.nodes <- dplyr::n_distinct(reloc_top_rssi$NodeId)
  
  # Determine the node with the strongest RSS value
  max.RSS <- reloc_top_rssi %>% slice_max(maxRSS_cs, n = 1, with_ties = FALSE)
  
  
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
    # extract location
    par.est <- cbind(coef(nls.test))
    
    # Estimated location of the point
    estimated.loc <- data.frame(
      reloc = reloc_id,
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
join_node_relocs <- node_relocs %>% select(reloc, TagId, Time.local, bin.2h, tag_type, tag_attachment_method,
                                         attachment_date, species, id, attachment_hour, Tag.start.time) %>%
  group_by(reloc) %>%
  slice(which.min(Time.local)) %>%
  distinct() %>%
  ungroup()

trilat_rssi_result <- left_join(trilat_rssi, join_node_relocs, by = "reloc")

write_csv(trilat_rssi_result, "trilat/full.csv")

## 4b. Second smoothing ------------------------------------------------------

# Cubic splines
cubic_spline_smooth <- do.call("rbind", lapply(unique(trilat_rssi_result$TagId), function(tag) {
  
  # select tag
    track <- trilat_rssi_result %>% filter(TagId == tag)
  
    bin_loop <- do.call("rbind", lapply(unique(track$bin.2h), function(bin){
      
      # select 2 h bin
      bin_df <- track %>% filter(bin.2h == bin)
    
    if (nrow(bin_df) >= 5) {
      
      sp_x <- smooth.spline(x = bin_df$Time.local, y = bin_df$x.est)
      sp_y <- smooth.spline(x = bin_df$Time.local, y = bin_df$y.est)
      
      data.frame(TagId = bin_df$TagId,
                 reloc = bin_df$reloc,
                 Time.local = bin_df$Time.local, 
                 x_cs = sp_x[["y"]], # smoothed value of x
                 y_cs = sp_y[["y"]]) # smoothed value of y
      
    } else {
      
      data.frame(TagId = bin_df$TagId,
                 reloc = bin_df$reloc,
                 Time.local = bin_df$Time.local,
                 x_cs = bin_df$x.est, # smoothed value of x
                 y_cs = bin_df$y.est) # smoothed value of y
      
    }
  }))
}))

## append results to dataframe
smoothed_tracks <- left_join(trilat_rssi_result, cubic_spline_smooth, by = c("TagId", "reloc", "Time.local"))

# save
write_csv(smoothed_tracks, "trilat/cs_smoothed_localizations.csv")


## 5. Signal det thru time ---------------------------------------------------

# tracks
track_time_df <- smoothed_tracks %>% mutate(radio_val = "Location estimate",
                                            duration = (as.numeric(Time.local - Tag.start.time)/86400)) %>%
  select(c(duration, radio_val, TagId))


time_df <- sel_df %>% mutate(radio_val = "Radio detection",
                             duration = (as.numeric(Time.local - Tag.start.time)/86400)) %>%
  select(c(duration, radio_val, TagId))

g_df <- bind_rows(track_time_df, time_df)

# give species info
df_sp_id <- data.frame(TagId = c("0752194C","1955342A","2D52614C",
                              "3366522D","3419331E","34522D78",
                              "4B071E66","4B4C4B2A","524C332D",
                              "4C4B1E2D"),
                    sp = c("Pterophanes_cyanopterus", "Pterophanes_cyanopterus", "Chalcostigma_heteropogon",
                           "Chalcostigma_heteropogon", "Chalcostigma_heteropogon", "Chalcostigma_heteropogon",
                           "Chalcostigma_heteropogon", "Pterophanes_cyanopterus", "Pterophanes_cyanopterus",
                           "Chalcostigma_heteropogon"))
g_df <- left_join(g_df, df_sp_id, by = "TagId")

g_df$radio_val <- factor(g_df$radio_val, levels = c("Location estimate", "Radio detection"))

# ordering from shortest to longest det duration by tagid
g_df$TagId <- factor(g_df$TagId, levels = c("3419331E", "2D52614C", "4B071E66", "3366522D", "4C4B1E2D", "34522D78",
                                            "524C332D","1955342A","4B4C4B2A","0752194C"))

(time_g <- ggplot() +
    geom_tile(data = g_df,
              aes(x = duration, y = radio_val, 
                  color = interaction(radio_val,sp)), 
               height = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c("#DB9D80", "#A0522D", "#00F5F5", "#008B8B")) +
    #scale_color_manual(values = c("#E59500", "#440154")) +
    facet_wrap(~TagId, scales = "free_x", nrow = 2, ncol = 6)+
    labs(y = "", x = "Days") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank()))

# total tracking and signal time
max_times <- g_df %>% group_by(TagId, radio_val) %>%
  summarise(max_days = max(duration))

## 6. Track cleaning ------------------------------------------------------

# import again if needed
#smoothed_tracks <- read_csv("trilat/cs_smoothed_localizations.csv")

# (tz_check <- ggplot() +
#     geom_density(data = smoothed_tracks, aes(x = hour(Time.local))))
# 
# smoothed_tracks$Time.local <- with_tz(smoothed_tracks$Time.local, tzone = "America/Bogota")

## time interval ---------------------------------------------------------------
# interval of estimated relocs, compare to initial radio signal interval
proj_crs <- st_crs("EPSG:3116") # Bogota Zone, already projected

smooth_sf <- st_as_sf(smoothed_tracks, coords = c("x_cs", "y_cs"),
                      crs = proj_crs)

distances_sf <- smooth_sf %>% 
  group_by(TagId) %>% # for each ind
  arrange(Time.local) %>% # order by time
  mutate(distance = st_distance(x = geometry,
                                y = lag(geometry),
                                by_element = TRUE)) %>%
  st_drop_geometry()%>%
  select(TagId, reloc, distance)

# make numeric
distances_sf$distance <- as.numeric(distances_sf$distance)

# replace NA values starting points
distances_sf[is.na(distances_sf$distance),]$distance <- 0

# join 
sm_tracks <- left_join(smoothed_tracks, distances_sf, by = c("TagId", "reloc"))

dt_tracks <- sm_tracks %>% 
  group_by(TagId) %>% # for each ind
  arrange(Time.local) %>% # order by time
  mutate(Time.dif = Time.local - lag(Time.local), # in seconds
         speed = distance/as.numeric(Time.dif), # m/s
         Day.start = as.POSIXct(strptime(attachment_date, format = "%m/%d/%Y", tz = "America/Bogota")),
         Date = format(Time.local, format="%m/%d/%Y"),
         Day.from.start = as.numeric(difftime(as.POSIXct(strptime(Date, format = "%m/%d/%Y", 
                                                                  tz = "America/Bogota")),
                                   Day.start, units = "days"))) %>%
  ungroup()

# NA time diff and speed are initial times
dt_tracks[is.na(dt_tracks$Time.dif),]$Time.dif <- as.duration(seconds(0))
dt_tracks[is.na(dt_tracks$speed),]$speed <- 0

unique(is.na(dt_tracks))

# mean dt by tag
dt_sum <- dt_tracks %>%
  group_by(TagId) %>%
  summarise(mean_dt = mean(Time.dif),
            med_dt = median(Time.dif),
            sum_days = max(Day.from.start))

quantile(dt_tracks$speed, probs = c(0.05, 0.5, 0.75, 0.9, 0.95)) # 95% under 24m/s, 90 % under 10 m/s, 75% under 2 m/s

# visualize speed distributions
(spped_g <- ggplot() +
    geom_density(data = dt_tracks, aes(x = speed)) +
    xlim(0, 24) +
    facet_wrap(~TagId, scales = "free"))

dt_tracks_filtered <- dt_tracks %>% filter(speed <= 25)


## 7. Autocorrelated kernels ------------------------------------------------------
# continuous-time movement models

## load necessary packages
library(tidyverse)
library(move)
library(ctmm)
library(sf)
library(ggspatial)
library(viridis)
library(cowplot)

# guide on using move objects:
# https://cran.r-project.org/web/packages/move/vignettes/move.html#import-non-movebank-data

# tutorials for ctmm:
# https://ctmm-initiative.github.io/ctmm/index.html
# https://ecoisilva.github.io/AKDE_minireview/code/AKDE_R-tutorial.html#step-3.-selecting-the-best-fit-movement-model-through-model-selection

## Do by tag
unique(dt_tracks_filtered$TagId)
id <- "1955342A"

# get data into move format
# first order according to individual and time stamp
ordered_tracks <- dt_tracks_filtered %>% filter(TagId == "1955342A") %>% arrange(Time.local)

move_df <- move(x = ordered_tracks$x_cs, y = ordered_tracks$y_cs,
                time = ordered_tracks$Time.local,
                proj = CRS("EPSG:3116"),
                data = ordered_tracks,
                animal = ordered_tracks$TagId)

# make ctmm object
tel <- as.telemetry(move_df)

# assign error
uere(tel) <- 97 # median error flying trials

# fit variogram
svf_1 <- variogram(tel)

# guess model estimates
guess_1 <- ctmm.guess(tel,
                      variogram = svf_1,
                      interactive = FALSE)

# fit model
fit_1 <- ctmm.select(tel,
                     CTMM = guess_1,
                     verbose = TRUE) # shows comparison of different models


summary(fit_1)

# calculate akde
akde <- akde(tel, fit_1)
summary(akde)$CI


# to spatial polygons dataframe for better graphing
akde_sp <- SpatialPolygonsDataFrame.UD(akde, level.UD = c(0.5, 0.75, 0.95))
akde_sf <- st_as_sf(akde_sp)
akde_sf$TagId <- id

# save
st_write(akde_sf, paste0("kernels/", id, ".shp"))

# clear up some space before going on to next one
rm(ordered_tracks)
rm(move_df)
rm(tel)
rm(svf_1)
rm(guess_1)
rm(fit_1)
rm(akde)
gc()

# read in shps
shp_files <- list.files("kernels/")
shp_files <- shp_files[grep("\\.shp$", shp_files)]

shp_list <- lapply(shp_files, function(file) {
  shp <- st_read(paste0("kernels/", file))
  return(shp)
})

shp_df <- bind_rows(shp_list)

shp_tidy <- shp_df %>%
  mutate(thresh = str_replace(name, "^.*?(50% low)", "50% low") %>%
           str_replace("^.*?(50% est)", "50% est") %>%
           str_replace("^.*?(50% high)", "50% high") %>%
           str_replace("^.*?(75% low)", "75% low") %>%
           str_replace("^.*?(75% est)", "75% est") %>%
           str_replace("^.*?(75% high)", "75% high") %>%
           str_replace("^.*?(95% low)", "95% low") %>%
           str_replace("^.*?(95% est)", "95% est") %>%
           str_replace("^.*?(95% high)", "95% high"))


dt_tracks_filtered$TagId <- factor(dt_tracks_filtered$TagId, 
                                   levels = c("2D52614C", "3366522D","3419331E","34522D78", "4B071E66", # Chalcostigma
                                              "0752194C","1955342A","4B4C4B2A","524C332D")) # Pterophanes
shp_tidy$TagId <- factor(shp_tidy$TagId, 
                         levels = c("2D52614C", "3366522D","3419331E","34522D78", "4B071E66", # Chalcostigma
                                    "0752194C","1955342A","4B4C4B2A","524C332D")) # Pterophanes

# option of taking away outlier
shp_tidy_f <- shp_tidy %>% filter(!TagId == "4B071E66")


# separate
est_50 <- shp_tidy_f %>% filter(thresh == "50% est") %>% mutate(thresh = "50%") # labels for pretty graph
lo_50 <- shp_tidy_f %>% filter(thresh == "50% low")
hi_50 <- shp_tidy_f %>% filter(thresh == "50% high")
est_75 <- shp_tidy_f %>% filter(thresh == "75% est") %>% mutate(thresh = "75%")
lo_75 <- shp_tidy_f %>% filter(thresh == "75% low")
hi_75 <- shp_tidy_f %>% filter(thresh == "75% high")
est_95 <- shp_tidy_f %>% filter(thresh == "95% est") %>% mutate(thresh = "90%")
lo_95 <- shp_tidy_f %>% filter(thresh == "95% low")
hi_95 <- shp_tidy_f %>% filter(thresh == "95% high")

# grid info
nodes <- read.csv("~/Mac_copy/CTT/new_own_approach/calibration_trial_tracks/Nodes.csv")
nodes_sel <- c("N01", "N02", "N03", "N04", "N08", "N44", "N45", "N46",
               "N43", "N40", "N23", "N30", "N37", "N36", "N35", "N34",
               "N33", "N32", "N31", "N24", "N17", "N13", "N09", "N05")

nodes_edge <- nodes %>% filter(Own_NodeId %in% nodes_sel)

# order sequentially
nodes_edge_ordered <- nodes_edge %>% arrange(match(Own_NodeId, nodes_sel))
first_point <- nodes_edge_ordered %>% slice_head() # add first node as last to close polygon
nodes_edge_closed <- bind_rows(nodes_edge_ordered, first_point)

# create line and polygon
nodes_line <- st_linestring(matrix(c(nodes_edge_closed$NodeUTMx, nodes_edge_closed$NodeUTMy), ncol = 2))
nodes_pol <- st_polygon(list(nodes_line))
proj_crs <- st_crs("EPSG:3116") # Bogota Zone, already projected
pol_sf <- st_sfc(nodes_pol, crs = proj_crs)


(hr_k_g <- ggplot() +
    geom_sf(data = est_95, aes(geometry=geometry, fill = thresh), color = NA, alpha = 0.7, show.legend = TRUE) +
    geom_sf(data = lo_95, aes(geometry=geometry), fill = NA, color = "#440154") +
    geom_sf(data = hi_95, aes(geometry=geometry), fill = NA, color = "#440154") +
    geom_sf(data = est_75, aes(geometry=geometry, fill = thresh), color = NA, alpha = 0.7, show.legend = TRUE) +
    geom_sf(data = lo_75, aes(geometry=geometry), fill = NA, color = "#21918c") +
    geom_sf(data = hi_75, aes(geometry=geometry), fill = NA, color = "#21918c") +
    geom_sf(data = est_50, aes(geometry=geometry, fill = thresh), color = NA, alpha = 0.7, show.legend = TRUE) +
    geom_sf(data = lo_50, aes(geometry=geometry), fill = NA, color = "#fde725") +
    geom_sf(data = hi_50, aes(geometry=geometry), fill = NA, color = "#fde725") +
    geom_sf(data = pol_sf, aes(geometry = geometry), color = NA, fill = c("#D4D4D4"), alpha = 0.7) +
    scale_fill_manual(values = c("#fde725", "#21918c", "#440154")) +
    facet_wrap(~TagId, nrow = 2) +
    scale_x_continuous(labels = c("","73.785°W", "","73.775°W", "")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.background=element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank()))

# calculate areas
st_area(pol_sf)
areas <- st_area(shp_tidy)
shp_tidy$area <- as.numeric(areas*0.000001) # in sq km

df_sp <- data.frame(TagId = c("0752194C","1955342A","2D52614C",
                              "3366522D","3419331E","34522D78",
                              "4B071E66","4B4C4B2A","524C332D"),
                    sp = c("Pterophanes_cyanopterus", "Pterophanes_cyanopterus", "Chalcostigma_heteropogon",
                           "Chalcostigma_heteropogon", "Chalcostigma_heteropogon", "Chalcostigma_heteropogon",
                           "Chalcostigma_heteropogon", "Pterophanes_cyanopterus", "Pterophanes_cyanopterus"))
k <- left_join(shp_tidy, df_sp, by = "TagId")

# arrange lo, est and hi estimates
k_est <- k %>% filter(thresh == "50% est" | thresh == "75% est" | thresh == "95% est") %>%
  select(TagId, sp, thresh, area) %>%
  st_drop_geometry() %>%
  mutate(est = area, 
         percent = case_when(thresh == "50% est" ~ "50%",
                             thresh == "75% est" ~ "75%",
                             thresh == "95% est" ~ "95%")) %>%
  select(-c(thresh, area))


k_low <- k %>% filter(thresh == "50% low" | thresh == "75% low" | thresh == "95% low") %>%
  select(TagId, sp, thresh, area) %>%
  st_drop_geometry() %>%
  mutate(lower_bound = area, 
         percent = case_when(thresh == "50% low" ~ "50%",
                             thresh == "75% low" ~ "75%",
                             thresh == "95% low" ~ "95%")) %>%
  select(-c(thresh, area))

k_high <- k %>% filter(thresh == "50% high" | thresh == "75% high" | thresh == "95% high") %>%
  select(TagId, sp, thresh, area) %>%
  st_drop_geometry() %>%
  mutate(upper_bound = area, 
         percent = case_when(thresh == "50% high" ~ "50%",
                             thresh == "75% high" ~ "75%",
                             thresh == "95% high" ~ "95%")) %>%
  select(-c(thresh, area))

k_tidy <- left_join(k_est, k_low, by = c("TagId", "sp", "percent"))
k_tidy <- left_join(k_tidy, k_high, by = c("TagId", "sp", "percent"))

k_tidy_f <- k_tidy %>% filter(!TagId == "4B071E66")
out_k <- k_tidy %>% filter(TagId == "4B071E66")

(bx_areas <- ggplot() +
    geom_pointrange(data = k_tidy, aes(x = percent, y = est, ymin = lower_bound, ymax= upper_bound, color = sp), 
                    show.legend = FALSE, position = position_jitterdodge(jitter.width = 0.1)) +
    geom_boxplot(data = k_tidy_f, aes(y = est, x = percent, color = sp), fill = NA, show.legend = FALSE) +
    scale_color_manual(values = c("#A0522D", "#008B8B")) +
    labs(x = "UD", y = "Area (km^2)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background=element_blank()))

mean_areas <- k_est %>% group_by(sp, percent) %>%
  summarise(mean_area = mean(est),
            sd_area = sd(est))

## map home range outlier
# separate
est_50 <- shp_tidy %>% filter(thresh == "50% est") %>% mutate(thresh = "50%") # labels for pretty graph
lo_50 <- shp_tidy %>% filter(thresh == "50% low")
hi_50 <- shp_tidy %>% filter(thresh == "50% high")
est_75 <- shp_tidy %>% filter(thresh == "75% est") %>% mutate(thresh = "75%")
lo_75 <- shp_tidy %>% filter(thresh == "75% low")
hi_75 <- shp_tidy %>% filter(thresh == "75% high")
est_95 <- shp_tidy %>% filter(thresh == "95% est") %>% mutate(thresh = "90%")
lo_95 <- shp_tidy %>% filter(thresh == "95% low")
hi_95 <- shp_tidy %>% filter(thresh == "95% high")

(hr_k_g_out <- ggplot() +
    geom_sf(data = est_95[est_95$TagId == "4B071E66",], aes(geometry=geometry, fill = thresh), color = NA, alpha = 0.7, show.legend = TRUE) +
    geom_sf(data = lo_95[lo_95$TagId == "4B071E66",], aes(geometry=geometry), fill = NA, color = "#440154") +
    geom_sf(data = hi_95[hi_95$TagId == "4B071E66",], aes(geometry=geometry), fill = NA, color = "#440154") +
    geom_sf(data = est_75[est_75$TagId == "4B071E66",], aes(geometry=geometry, fill = thresh), color = NA, alpha = 0.7, show.legend = TRUE) +
    geom_sf(data = lo_75[lo_75$TagId == "4B071E66",], aes(geometry=geometry), fill = NA, color = "#21918c") +
    geom_sf(data = hi_75[hi_75$TagId == "4B071E66",], aes(geometry=geometry), fill = NA, color = "#21918c") +
    geom_sf(data = est_50[est_50$TagId == "4B071E66",], aes(geometry=geometry, fill = thresh), color = NA, alpha = 0.7, show.legend = TRUE) +
    geom_sf(data = lo_50[lo_50$TagId == "4B071E66",], aes(geometry=geometry), fill = NA, color = "#fde725") +
    geom_sf(data = hi_50[hi_50$TagId == "4B071E66",], aes(geometry=geometry), fill = NA, color = "#fde725") +
    geom_sf(data = pol_sf, aes(geometry = geometry), color = NA, fill = c("#D4D4D4"), alpha = 0.7) +
    scale_fill_manual(values = c("#fde725", "#21918c", "#440154")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.background=element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank()))
