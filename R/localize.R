## Update localizaitons using multilateration method
library(foreach)
library(dplyr)

## Localize function
source("./R/localize_fn.R")

cores <- parallel::detectCores()
cl <- parallel::makeForkCluster(cores-11, outfile = "")
doParallel::registerDoParallel(cl)

## Prepared files
files = list.files("./outputs/prepared/")

foreach(f = files,
        .packages=c("tidyverse","geosphere"),
        .verbose = TRUE) %dopar%
  { localize_fn(prepared_file =  f,
                reps = 100,
                model_scale = "exp")
  }

