## Update localizaitons using multilateration method
library(foreach)
library(dplyr)

## Localize function
source("./R/localize_fn.R")

cl <- parallel::makeForkCluster(3, outfile = "")
doParallel::registerDoParallel(cl)

## Prepared files
files = list.files("./outputs/prepared/")

## Model path
rssi_dist_model_file = "./data/RSSI_log_dist_model_zebby.RDS"

foreach(f = files,
        .packages=c("tidyverse","geosphere"),
        .verbose = TRUE) %dopar%
  { localize_fn(prepared_file =  f,
                rssi_dist_model_file = rssi_dist_model_file,
                reps = 100)
  }

