## CALIBRATION OF CTT DATA ##
## Compare effects of different factors
# on calculated error

# load required packages
library(tidyverse)
library(nlme)
library(DHARMa)
library(MuMIn)
library(lubridate)
library(cowplot)
library(sf)
library(raster)
library(ggspatial)
library(viridis)

## 1. Calculate summary error var -----------------------------------

df_all <- read.csv("df_all_models_error.csv")

# by trials
sum.stats_trials <- df_all %>%
  group_by(type_signal, loc_method, param,
           smoothing, second_smoothing, test_height, unique_trial_id) %>%
  summarise(min.diff_2.5 = quantile(error_m, probs = 0.025),
            max.diff_97.5 = quantile(error_m, probs = 0.975),
            sd.diff = sd(error_m),
            med.diff = median(error_m),
            avg.diff = mean(error_m),
            h.avg.diff = (1/mean(1/error_m)),
            n.relocs = n()) %>%
  ungroup()

mean(sum.stats_trials$sd.diff, na.rm = TRUE)
sd(sum.stats_trials$sd.diff, na.rm = TRUE)

# remove trials with less than one estimated reloc
sum.stats_trials_sel <- sum.stats_trials %>% filter(n.relocs > 1)

sum.stats_gen <- sum.stats_trials_sel %>%
  group_by(type_signal, loc_method, param,
           smoothing, second_smoothing, test_height) %>%
  summarise(mean_min.diff_2.5 = mean(min.diff_2.5),
            mean_max.diff_97.5 = mean(max.diff_97.5),
            mean_sd.diff = mean(sd.diff),
            mean_med.diff = mean(med.diff),
            mean_avg.diff = mean(avg.diff),
            mean_h.avg.diff = mean(h.avg.diff),
            mean_n.relocs = mean(n.relocs),
            total_n.relocs = sum(n.relocs),
            n.trials = n()) %>%
  ungroup()


## 2. LMM: process workflow -------------------------------

## Tidy data ---------------------------------------------------------
# remove tests with only one relocation 
df_count_single_relocs <- df_all %>% group_by(unique_trial_id, test_height, loc_method,
                                          param, type_signal, smoothing, second_smoothing) %>%
  summarise(n_test = n())

single_relocs <- df_count_single_relocs %>% filter(n_test == 1)
single_trial_id <- unique(single_relocs$unique_trial_id) # 4 trials

df_model <- df_all %>% filter(!unique_trial_id %in% single_trial_id)

df_check <- df_model %>% group_by(unique_trial_id, test_height, loc_method,
                                          param, type_signal, smoothing, second_smoothing) %>%
  summarise(n_test = n())
range(df_check$n_test)

# all possible combinations for order and ACF calc later
unique_trials <- unique(df_model$unique_trial_id)
unique_height <- unique(df_model$test_height)
unique_loc <- unique(df_model$loc_method)
unique_param <- unique(df_model$param)
unique_signal <- unique(df_model$type_signal)
unique_smooth <- unique(df_model$smoothing)
unique_second_smooth <- unique(df_model$second_smoothing)

combinations <- expand.grid(
  loc_method = unique_loc,
  param = unique_param,
  type_signal = unique_signal,
  smoothing = unique_smooth,
  second_smoothing = unique_second_smooth
)

# arrange df by time
comb_arrange <- combinations %>% mutate(method_id = seq(1, nrow(combinations)))
comb_df <- left_join(df_model, comb_arrange, by = c("loc_method", "param","type_signal", 
                                                    "smoothing", "second_smoothing"))
df_ordered <- comb_df %>% arrange(method_id, test_local_time)

# force baseline factor levels for model comparison

# loc_method = multilateration
# param = general
# type_signal = avgRSS
# smoothing = none

# type signal average already selected since it starts with "a"
# and second smoothing since "n" is before ss

df_ordered[df_ordered$loc_method == "Multilateration",]$loc_method <- "a_Multilateration"
df_ordered[df_ordered$param == "General",]$param<- "a_General"
df_ordered[df_ordered$smoothing == "None",]$smoothing<- "a_None"
df_ordered[df_ordered$second_smoothing == "ss_None",]$second_smoothing<- "a_None"

# check distribution of error
hist(df_model$error_m)
hist(sqrt(df_model$error_m))

## Fit model ---------------------------------------------------------

full_lm <- lme(sqrt(error_m) ~
                  loc_method +
                  param +
                  type_signal +
                  smoothing +
                  second_smoothing,
                  random = ~ 1|unique_trial_id/reloc_id,
                data = df_ordered)


summary(full_lm)

resid_acf <- ACF(full_lm, resType = "normalized")
resid_acf_plot <- plot(resid_acf)

# residuals and model fit
resid <- residuals(full_lm, type = "normalized")
df_ordered$resid <- resid

# calculate residual in next time step
df_model_sorted <- df_ordered %>% group_by(method_id) %>%
  arrange(test_local_time) %>%
  mutate(next_resid = lag(resid))


# visualize
(autocor_plot <- ggplot() +
    geom_point(data = df_model_sorted, aes(x = resid, y = next_resid), size = 1, alpha = 0.3, color = "#454545") +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Residual (t)", y = "Residual (t+1)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))



## fit model without random effects to compare temporal autocorrelation
no_re_lm <- gls(sqrt(error_m) ~ 
                 loc_method +
                 param +
                 type_signal +
                 smoothing +
                 second_smoothing,
               data = df_ordered)

resid_re_acf <- ACF(no_re_lm, resType = "normalized")
resid_re_plot <- plot(resid_re_acf)

# residuals and model fit
resid_nore <- residuals(no_re_lm, type = "normalized")
df_ordered$resid_nore <- resid_nore

# calculate residual in next time step
df_model_sorted_nore <- df_ordered %>% group_by(method_id) %>%
  arrange(test_local_time) %>%
  mutate(next_resid_nore = lag(resid_nore))

# visualize
(autocor_re_plot <- ggplot() +
    geom_point(data = df_model_sorted_nore, aes(x = resid_nore, y = next_resid_nore), size = 1, alpha = 0.3, color = "#454545") +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Residual (t)", y = "Residual (t+1)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

# plot for comparison
plot_grid(resid_re_plot, autocor_re_plot,
          resid_acf_plot, autocor_plot, 
          nrow = 2, labels = "AUTO")

## Model selection ---------------------------------------------------------
options(na.action = 'na.fail')

d_df <- dredge(full_lm, evaluate = TRUE, rank = "AICc")


## Compare fixed effects ---------------------------------------------------------

# get confidence intervals
ci <- intervals(full_lm)
ci_df <- as.data.frame(ci$fixed)
ci_df$fixed_effect <- rownames(ci_df)

# remove intercept
ci_df <- ci_df %>% filter(!fixed_effect == "(Intercept)")

# rename
ci_df <- ci_df %>% mutate(fixed_effect = case_when(fixed_effect == "loc_methodK-means" ~ "K-means",
                                                   fixed_effect == "loc_methodNearest nodes" ~ "Nearest nodes",
                                                   fixed_effect == "loc_methodStrongest signal" ~ "Strongest signal",
                                                   fixed_effect == "paramBy node" ~ "By node",
                                                   fixed_effect == "paramBy tag" ~ "By tag",
                                                   fixed_effect == "type_signalmaxRSS" ~ "maxRSS",
                                                   fixed_effect == "smoothingCubic splines" ~ "Signal cubic splines",
                                                   fixed_effect == "smoothingKalman" ~ "Signal Kalman",
                                                   fixed_effect == "second_smoothingss_Cubic" ~ "Track cubic splines",
                                                   fixed_effect == "second_smoothingss_Kalman" ~ "Track Kalman"))
# add decision groups
ci_df$group <- c("NLS multilateration", "NLS multilateration", "NLS multilateration", 
                 "Model parameters", "Model parameters", "Signal type",
                 "Signal smoothing", "Signal smoothing",
                 "Track smoothing", "Track smoothing")

# add baselines
#ci_df <- ci_df[,4:8] # remove unnecessary cols
baselines <- data.frame(fixed_effect = c("All nodes", "General", "avgRSS", "No signal smooth", "No track smooth"),
                        group = c("NLS multilateration", "Model parameters", "Signal type", "Signal smoothing", "Track smoothing"),
                        lower = 0,
                        est. = 0,
                        upper = 0)

g_df <- bind_rows(ci_df, baselines)

# add increase/neutral/decrease to color graph
g_df$change_effect <- c("d", "i", "d", "i", "i", "d", "d", "d", "d", "d", "n", "n", "n", "n", "n")

# organize order
g_df$group <- factor(g_df$group, levels = c("Signal type", "Signal smoothing", "Model parameters", "NLS multilateration", "Track smoothing"))

g_df$fixed_effect <- factor(g_df$fixed_effect , levels = c("maxRSS","avgRSS",
                                         "Signal Kalman", "Signal cubic splines", "No signal smooth",
                                         "By tag", "By node", "General",
                                         "K-means", "Nearest nodes", "Strongest signal", "All nodes",
                                         "Track Kalman", "Track cubic splines", "No track smooth"))

# plot
(g_error_lm <- ggplot(data = g_df) +
    geom_vline(xintercept = 0, linetype="dashed", color = "#FFB90F") +
    geom_point(aes(y = fixed_effect, x = est., color = change_effect), 
               show.legend = FALSE) +
    geom_linerange(aes(y = fixed_effect, xmin = lower, xmax= upper, color = change_effect), 
                   show.legend = FALSE) +
    scale_color_manual(values = c("#436EEE", "#FF3E96", "#FFB90F")) +
    facet_wrap(~group, scales = "free_y", ncol = 1) +
    labs(y = "", x = "Estimated effect") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()))


## 3. LMM: environmental var ---------------------------------------------------

## Tidy data ---------------------------------------------------------

# model on selected method
model_method_sel <- df_ordered %>% filter(type_signal == "maxRSS" & loc_method == "Strongest signal" &
                                  param == "By node" & smoothing == "Cubic splines" &
                                  second_smoothing == "ss_Cubic")


# assign vegetation covers
veg_layers <- st_read("vegetation_layers/grid_map/grid_50buffer_1mres.shp")

proj_crs <- st_crs("EPSG:3116") # Bogota Zone
points <- st_as_sf(model_method_sel, coords = c("TestUTMx", "TestUTMy"),
                   crs = proj_crs)

points_buffer <- points %>% st_buffer(100)

points_veg <- st_join(points_buffer, veg_layers)

f <- veg_layers %>% filter(value == "Forest" | value == "Bamboo")
hp <- veg_layers %>% filter(value == "High_paramo")
op <- veg_layers %>% filter(value == "Open_paramo")
r <- veg_layers %>% filter(value == "Road")

veg_p <- do.call("rbind", lapply(unique(points_buffer$reloc_id), function(reloc) {
  
  p <- points_buffer %>% filter(reloc_id == reloc)
  
  f_c <- st_crop(f, p)
  f_a <- st_area(f_c)
  
  hp_c <- st_crop(hp, p)
  hp_a <- st_area(hp_c)
  
  op_c <- st_crop(op, p)
  op_a <- st_area(op_c)
  
  r_c <- st_crop(r, p)
  r_a <- st_area(r_c)
  
  
  data.frame(reloc_id = p$reloc_id,
             f_area = sum(f_a),
             hp_area = sum(hp_a),
             op_area = sum(op_a),
             r_area = sum(r_a))
  
}))

# this took a long time to run,
# save just in case
write_csv(veg_p, "veg_p_temporal.csv")


df_model_env <- left_join(model_method_sel, veg_p, by = "reloc_id")

# calculate percentages within circle
df_model_env <- df_model_env %>% mutate(f_percent = as.numeric(f_area)/31401.57,
                                        hp_percent = as.numeric(hp_area)/31401.57,
                                        op_percent = as.numeric(op_area)/31401.57,
                                        r_percent = as.numeric(r_area)/31401.57) %>%
  mutate(f_percent = case_when(f_percent > 1 ~ 1,
                               f_percent < 1 ~ f_percent),
         hp_percent = case_when(hp_percent > 1 ~ 1,
                               hp_percent < 1 ~ hp_percent),
         op_percent = case_when(op_percent > 1 ~ 1,
                               op_percent < 1 ~ op_percent),
         r_percent = case_when(r_percent > 1 ~ 1,
                               r_percent < 1 ~ r_percent))

# road is not a dominant cover for any point
max(df_model_env$r_percent) # max 0.09

# for each point, select dominating veg cover
df_model_env_long <- df_model_env %>% pivot_longer(cols = c("f_percent", "hp_percent", "op_percent"),
                                                   names_to = "veg_type", values_to = "percent")

dominant_veg <- do.call("rbind", lapply(unique(df_model_env_long$reloc_id), function(reloc) {
  
  rel <- df_model_env_long %>% filter(reloc_id == reloc)
  
  dominant <- rel %>% slice_max(percent, with_ties = FALSE)
  
  data.frame(reloc_id = rel$reloc_id,
             dominant_veg_type = dominant$veg_type,
             dominant_veg_type_percent = dominant$percent)
  
}))

dominant_veg_j <- dominant_veg %>% distinct()

df_model_env <- left_join(df_model_env, dominant_veg_j, by = "reloc_id")


# assign topography
# downloaded using following link of region DSM_10_N04_00_W074_00
# https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-30-DGED__2022_1/Copernicus_DSM_10_N04_00_W074_00.tar
dem <- raster("Copernicus_DEM/Copernicus_DSM_10_N04_00_W074_00/Copernicus_DSM_10_N04_00_W074_00/DEM/Copernicus_DSM_10_N04_00_W074_00_DEM.tif")

# crop to area of interest
# spatial points back to latlon match dem proj

latlon_crs <- st_crs("EPSG:4326")
bbox <- points %>% st_buffer(200) %>%
  st_transform(crs = latlon_crs) %>% st_bbox()

area <- extent(c(-73.774941, -73.761176, 4.521257, 4.533782))
dem_c <- crop(dem, area)

# convert relocs to spatial points
# back to lat lon data
points_latlon <- points %>%
  st_transform(crs = latlon_crs)

points_100buffer_latlon <- points %>% st_buffer(100) %>%
  st_transform(crs = latlon_crs)

# Extract values of raster at points
ground_elevation <- extract(dem_c, points_latlon)
sd_ground_elevation <- extract(dem_c, points_100buffer_latlon, 
                               fun = sd)

df_model_env$ground_elevation <- ground_elevation
df_model_env$ruggedness <- sd_ground_elevation[,1]
df_model_env$flight_height_ground <- df_model_env$elevation - df_model_env$ground_elevation

# tidy up in or out of grid
df_model_env$in_or_out <- as.factor(df_model_env$in_or_out)

# save for later use
write_csv(df_model_env, "df_model_spatial.csv")

# df_model_env <- read_csv("df_model_spatial.csv")


## Vis distributions and correlations ------------------------------------------

# correlation
library(GGally)
quant_env <- df_model_env %>% select(dist_edge, ground_elevation, ruggedness, flight_height_ground)
ggpairs(quant_env)

# any outliers?
quantile(df_model_env$error_m, probs = c(0.025, 0.5, 0.85, 0.9, 0.95, 0.975))

df_model_env <- df_model_env %>% filter(error_m < 410 & error_m > 15)


# histograms
hist(df_model_env$error_m)
hist(sqrt(df_model_env$error_m))

hist(df_model_env$dist_edge)
hist(sqrt(df_model_env$dist_edge))

hist(df_model_env$ruggedness)
hist(sqrt(df_model_env$ruggedness))

hist(df_model_env$flight_height_ground)
hist(sqrt(df_model_env$flight_height_ground))

# scale and transform
df_model_env$dist_edge_s <- scale(sqrt(df_model_env$dist_edge))

df_model_env$ruggedness_s <- scale(sqrt(df_model_env$ruggedness))

df_model_env$flight_height_s <- scale(sqrt(df_model_env$flight_height_ground))

# check no nas
unique(is.na(df_model_env))


# vis possible relations
(g1 <- ggplot() +
    geom_point(data = df_model_env, aes(x = dist_edge, y = error_m)))
(g2 <- ggplot() +
    geom_point(data = df_model_env, aes(x = ruggedness, y = error_m)))
(g3 <- ggplot() +
    geom_point(data = df_model_env, aes(x = flight_height_ground, y = error_m)))
(g4 <- ggplot() +
    geom_boxplot(data = df_model_env, aes(x = dominant_veg_type, y = error_m)) +
    geom_point(data = df_model_env, aes(x = dominant_veg_type, y = error_m)))
(g5 <- ggplot() +
    geom_boxplot(data = df_model_env, aes(x = in_or_out, y = error_m)) +
    geom_point(data = df_model_env, aes(x = in_or_out, y = error_m)))

## Fit model ---------------------------------------------------------

df_model_env <- df_model_env %>% mutate(dominant_veg_type = 
                                          case_when(dominant_veg_type == "f_percent" ~ "Forest",
                                                    dominant_veg_type == "hp_percent" ~ "High_paramo",
                                                    dominant_veg_type == "op_percent" ~ "a_Open_paramo"))


# fit model of spatial features
spatial_lm <- lme(sqrt(error_m) ~
                    dist_edge_s +
                    in_or_out +
                    ruggedness_s +
                    flight_height_s +
                    # veg covers +
                    dominant_veg_type +
                    # interactions
                    dist_edge_s*dominant_veg_type +
                    in_or_out*dominant_veg_type +
                    ruggedness_s*dominant_veg_type +
                    flight_height_s*dominant_veg_type,
                  random = ~ 1|unique_trial_id,
                  data = df_model_env)

summary(spatial_lm)

# check acf
resid_spatial_acf <- ACF(spatial_lm, resType = "normalized")
resid_acf_spatial_plot <- plot(resid_spatial_acf)


spatial_lm_nore <- gls(sqrt(error_m) ~
                         dist_edge_s +
                         in_or_out +
                         ruggedness_s +
                         flight_height_s +
                         # veg covers +
                         dominant_veg_type +
                         # interactions
                         dist_edge_s*dominant_veg_type +
                         in_or_out*dominant_veg_type +
                         ruggedness_s*dominant_veg_type +
                         flight_height_s*dominant_veg_type,
                       data = df_model_env)

# check acf
resid_spatial_acf_nore <- ACF(spatial_lm_nore, resType = "normalized")
resid_acf_spatial_nore_plot <- plot(resid_spatial_acf_nore)

## Thin data to reduce ACF -----------------------------------------------------

df_env_thin <- do.call("rbind", lapply(unique(df_model_env$unique_trial_id), function(trial) {
  
  trial <- df_model_env %>% filter(unique_trial_id == trial)
  
  trial_thinned <- trial[seq(1, nrow(trial), 2), ]
  
}))

spatial_lm_thinned <- lme(sqrt(error_m) ~
                            dist_edge_s +
                            in_or_out +
                            ruggedness_s +
                            flight_height_s +
                            # veg covers +
                            dominant_veg_type +
                            # interactions
                            dist_edge_s*dominant_veg_type +
                            in_or_out*dominant_veg_type +
                            ruggedness_s*dominant_veg_type +
                            flight_height_s*dominant_veg_type,
                          random = ~ 1|unique_trial_id,
                          data = df_env_thin)

summary(spatial_lm_thinned)

# check acf
resid_spatial_acf_t <- ACF(spatial_lm_thinned, resType = "normalized")
resid_acf_spatial_plot_t <- plot(resid_spatial_acf_t)


## calculate residuals in next time step for all models
resid_env <- residuals(spatial_lm, type = "normalized")
df_model_env$resid <- resid_env

resid_env_nore <- residuals(spatial_lm_nore, type = "normalized")
df_model_env$resid_nore <- resid_env_nore

resid_thin <- residuals(spatial_lm_thinned, type = "normalized")
df_env_thin$resid <- resid_thin

# calculate residual in next time step
df_env_r <- df_model_env %>% group_by(unique_trial_id) %>%
  mutate(next_resid = lag(resid),
         next_resid_nore = lag(resid_nore))

df_env_thin_r <- df_env_thin %>% group_by(unique_trial_id) %>%
  mutate(next_resid = lag(resid))

# visualize
(plot_env_r <- ggplot() +
    geom_point(data = df_env_r, aes(x = resid, y = next_resid), size = 1, alpha = 0.3, color = "#454545") +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Residual (t)", y = "Residual (t+1)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

(plot_env_r_nore <- ggplot() +
    geom_point(data = df_env_r, aes(x = resid_nore, y = next_resid_nore), size = 1, alpha = 0.3, color = "#454545") +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Residual (t)", y = "Residual (t+1)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

(plot_env_r_thin <- ggplot() +
    geom_point(data = df_env_thin_r, aes(x = resid, y = next_resid), size = 1, alpha = 0.3, color = "#454545") +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Residual (t)", y = "Residual (t+1)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))


plot_grid(resid_acf_spatial_nore_plot, resid_acf_spatial_plot, resid_acf_spatial_plot_t,
          plot_env_r, plot_env_r_nore, plot_env_r_thin,
          nrow = 2, labels = c("A", "B", "C"))

## Model selection ---------------------------------------------------------
options(na.action = 'na.fail')

d_spatial <- dredge(spatial_lm_thinned, evaluate = TRUE, rank = "AICc")

# save
write_csv(as.data.frame(d_spatial),"spatial_influence_error_model_rank.csv")

# average first models 
sel <- d_spatial %>% filter(weight > 0.2) # 3 models, 77% AICc weights, delta from 0.99 to 2.73 (more than doubled)
sum(sel$weight)

m_avg_list <- get.models(d_spatial, subset = weight > 0.2)

# full average of model assumes coefficient = 0 where variable not present in model
summary(model.avg(m_avg_list))
confint(model.avg(m_avg_list))


## Compare fixed effects ---------------------------------------------------------

confint_env <- confint(model.avg(m_avg_list))
coef_env <- summary(model.avg(m_avg_list))[["coefficients"]]
coef_env <- coef_env[1,] #full averaged coef

param_env <- as.data.frame(confint_env)
param_env$estimate <- as.vector(coef_env)

# tidy names for fixed effects
param_env$fixed_effect <- c("Intercept",
                            "Edge distance",
                            "Dense veg.",
                            "Paramo",
                            "Flight height",
                            "Inside grid",
                            "Ruggedness",
                            "Dense veg. : Flight height",
                            "Paramo : Flight height",
                            "Dense veg. : Inside grid",
                            "Paramo : Inside grid",
                            "Dense veg. : Ruggedness",
                            "Paramo : Ruggedness",
                            "Dense veg. : Edge distance",
                            "Paramo : Edge distance")

# remove intercept
param_env <- param_env %>% filter(!fixed_effect == "Intercept")

# add change effects
param_env$change_effect <- c("n", "d", "d", "d", "d", "i", "n", "i", "i", "i", "n", "d", "n", "n")

# rename var
param_env <- param_env %>% rename(lower = "2.5 %", upper = "97.5 %")

# order
param_env$fixed_effect <- factor(param_env$fixed_effect, levels = c("Paramo : Flight height",
                                                                    "Dense veg. : Flight height",
                                                                    "Paramo : Ruggedness",
                                                                    "Dense veg. : Ruggedness",
                                                                    "Paramo : Inside grid",
                                                                    "Dense veg. : Inside grid",
                                                                    "Paramo : Edge distance",
                                                                    "Dense veg. : Edge distance",
                                                                    "Dense veg.",
                                                                    "Paramo",
                                                                    "Flight height",
                                                                    "Ruggedness",
                                                                    "Inside grid",
                                                                    "Edge distance"))

param_env$group <- "Spatial features"

# plot
(g_error_lm_env <- ggplot(data = param_env) +
   geom_vline(xintercept = 0, linetype="dashed", color = "#FFB90F") +
   geom_point(aes(y = fixed_effect, x = estimate, color = change_effect), 
              show.legend = FALSE) +
   geom_linerange(aes(y = fixed_effect, xmin = lower, xmax= upper, color = change_effect), 
                  show.legend = FALSE) +
   scale_color_manual(values = c("#436EEE", "#FF3E96", "#FFB90F")) +
   facet_wrap(~group) +
   labs(y = "", x = "Estimated effect") +
   theme_bw() +
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank()))


df_bx_graph <- df_model_env %>% mutate(in_or_out = case_when(in_or_out == 0 ~ "Outside",
                                              in_or_out == 1 ~ "Inside")) %>%
  mutate(dominant_veg_type = case_when(dominant_veg_type == "a_Open_paramo" ~ "Open",
                                       dominant_veg_type == "High_paramo" ~ "Paramo",
                                       dominant_veg_type == "Forest" ~ "Dense veg."))
df_bx_graph$dominant_veg_type <- factor(df_bx_graph$dominant_veg_type, 
                                        levels = c("Open", "Paramo", "Dense veg."))

(g_bx <- ggplot() +
    geom_boxplot(data = df_bx_graph, aes(x = in_or_out, y = error_m, fill = dominant_veg_type)) +
    scale_fill_manual(values = c("#FFF68F", "#D7E8BA", "#698B22")) +
    labs(y = "Error (m)", x = "Grid location", fill = "") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()))


## 4. Error quantile figure ----------------------------------------------------

# read simulated tracks
sim_df <- read_csv("ctmm_simulated_tracks_with_error.csv")

med_sim_df_error <- sim_df %>% group_by(test_height, unique_trial_id, reloc_within_trial) %>%
  summarise(med_error = median(error_m))

# calculate quantiles
q_df <- df_model_env %>% group_by(unique_trial_id) %>% # only selected method
  mutate(reloc_rank = rank(error_m),
         percent_rank = reloc_rank/n())

q_df_sim <- med_sim_df_error %>% group_by(unique_trial_id) %>% 
  mutate(reloc_rank = rank(med_error),
         percent_rank = reloc_rank/n())


# order facets
q_df$test_height <- factor(q_df$test_height, levels = c("drone", "mid_level", "ground_level"))
q_df_sim$test_height <- factor(q_df_sim$test_height, levels = c("drone", "mid_level", "ground_level"))

# labels for facets
test_labels <- c("High flight", "Low flight", "Ground")
names(test_labels) <- c("drone", "mid_level", "ground_level")


(q_g<- ggplot() +
    geom_line(data = q_df, aes(x = error_m, y = percent_rank, 
                                     group = unique_trial_id),
              color = "#454545") +
    geom_point(data = q_df, aes(x = error_m, y = percent_rank,
                                      color = ruggedness),
               alpha = 0.5) +
    geom_line(data = q_df_sim, aes(x = med_error, y = percent_rank, 
                               group = unique_trial_id),
              color = "#454545") +
    facet_wrap(~test_height, 
               labeller = labeller(test_height = test_labels),
               ncol = 1) +
    scale_color_viridis(option = "inferno") +
    labs(x = "Error (m)", y = "Quantile", color = "Ruggedness") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.title = element_text(size = 9)))


## 5. Difference with simulated tracks ----------------------------------------------------

# median error by trial for random tracks
med_sim_df_error <- sim_df %>% group_by(test_height, unique_trial_id) %>%
  summarise(sim_med_error = median(error_m)) %>%
  arrange(unique_trial_id)

# tidy up error info on all methods
unique_combinations <- expand.grid(
  loc_method = unique(df_all$loc_method),
  param = unique(df_all$param),
  type_signal = unique(df_all$type_signal),
  smoothing = unique(df_all$smoothing),
  second_smoothing = unique(df_all$second_smoothing)
)

# arrange df by time
un_comb <- unique_combinations %>% mutate(method_id = seq(1, nrow(unique_combinations)))
un_comb_df <- left_join(df_all, un_comb, by = c("loc_method", "param","type_signal", 
                                                    "smoothing", "second_smoothing"))
comb_vector <- unique(un_comb$method_id)

# loop over methods
sim_wilcox <- do.call("rbind", lapply(unique(un_comb$method_id), function(comb_id) {
  comb_error <- un_comb_df %>% filter(method_id == comb_id) %>%
    group_by(unique_trial_id) %>%
    summarise(med_error = median(error_m)) # median error by trial
  
  error_united <- left_join(comb_error, med_sim_df_error, by = "unique_trial_id")
  
  test <- wilcox.test(error_united$med_error, error_united$sim_med_error, paired = TRUE)
  
  data.frame(method_id = comb_id,
             test_statistic = test$statistic,
             p_val = test$p.value)
}))

wilcox_results_sum <- left_join(un_comb, sim_wilcox, by = "method_id")

range(wilcox_results_sum$p_val)