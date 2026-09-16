
source("00_moments.R")        # g_lognormal
source("00_helper_func.R")    # in_hull, etel_manual, etel_fit_manual, safe_fit


library(readr)
library(gmm)
library(tidyverse)
 

mu     <- 0
sigma2 <- 0.25
tet0   <- c(mu, sigma2)
 
 
 
all_data <- read_rds("outputs/02_all_data_sims.rds")
gel_results <- all_data |>
  group_by(pct, delta, sim) |>
  nest() |>
  ungroup() |>
  mutate(
    y = map(data, "y"),
    loglike_true = map_dbl(y, \(yy) etel_manual(tet0,yy)$loglike), #loglike at the true
    fit_ETEL = map(y, \(yy) safe_fit(tet0, yy)),
    status = if_else(map_lgl(fit_ETEL, is.null), "fit failed", "ok"),
    coef = map(fit_ETEL, "coef"),
    loglike = map_dbl(fit_ETEL, \(f) if (is.null(f)) NA_real_ else f$loglike), # loglike evaluated starting from the truth
    conv = map_int(fit_ETEL, \(f) if (is.null(f)) NA_integer_ else f$conv)
  )





# gel_results <- all_data |>
#   group_by(pct, delta, sim) |>
#   nest() |>
#   ungroup() |>
#   mutate(
#     y            = map(data, "y"),
#     loglike_true = map_dbl(y, \(yy) etel_manual(tet0, yy)$loglike),   # surface at the truth
#     fit_ETEL     = map(y, \(yy) safe_fit(tet0, yy)),
#     status       = if_else(map_lgl(fit_ETEL, is.null), "fit failed", "ok"),
#     coef         = map(fit_ETEL, "coef"),
#     loglike      = map_dbl(fit_ETEL, \(f) if (is.null(f)) NA_real_ else f$loglike),
#     conv         = map_int(fit_ETEL, \(f) if (is.null(f)) NA_integer_ else f$conv)
#   )
 
count(gel_results, pct, delta, status, conv)
 
write_rds(gel_results, "outputs/04_gel_results.rds", compress = "gz")

