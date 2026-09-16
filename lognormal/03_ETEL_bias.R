source("00_moments.R")
source("00_helper_func.R")



all_data <- read_rds("outputs/02_all_data_sims.rds")


# testing
# test_df <- all_data |>
#   filter(pct == 10, delta == 4, sim == 5)
# mu <- 0
# sigma2 <- 0.25
# biased_start(test_df$y)

mu     <- 0
sigma2 <- 0.25
frac   <- 0.7        
dir    <- c(1, 1) 

gel_results_bias <- all_data |>
  group_by(pct, delta, sim) |>
  nest() |>
  ungroup() |>
  mutate(
    y        = map(data, "y"),
    start    = map(y, \(yy) biased_start(yy, frac = frac, dir = dir)),
    tet0     = map(start, "tet0"),
    t_hull   = map_dbl(start, "t_hull"), # this is the distance from the tru parameter to the edge of the ETEL feasible region
    fit_ETEL = map2(tet0, y, safe_fit),
    status   = if_else(map_lgl(fit_ETEL, is.null), "fit failed", "ok"),
    coef     = map(fit_ETEL, "coef"),
    loglike  = map_dbl(fit_ETEL, \(f) if (is.null(f)) NA_real_ else f$loglike),
    conv     = map_int(fit_ETEL, \(f) if (is.null(f)) NA_integer_ else f$conv),
    frac     = frac
  ) |>
  select(-start)



count(gel_results_bias, pct, delta, status, conv)
 
write_rds(gel_results_bias, "outputs/04_gel_results_bias.rds", compress = "gz")



# # try this

# source("00_moments.R")   # g_lognormal
 
# library(readr)
# library(gmm)
# library(tidyverse)
 
 
# # -------------------- settings ----------------
 
# mu     <- 0
# sigma2 <- 0.25
# frac   <- 0.7          # fraction of the feasible distance from truth to hull boundary
# dir    <- c(1, 1)      # direction of the bias in (mu, sigma2)
 
 
# # -------------------- helpers ----------------
 
# in_hull <- function(tet, y) {
#   G  <- g_lognormal(tet, y)
#   lp <- lpSolve::lp(
#     direction = "min", objective.in = rep(0, nrow(G)),
#     const.mat = rbind(t(G), rep(1, nrow(G))),
#     const.dir = rep("=", ncol(G) + 1),
#     const.rhs = c(rep(0, ncol(G)), 1)
#   )
#   lp$status == 0
# }
 
# etel_manual <- function(tet, y) {
#   G         <- g_lognormal(tet, y)
#   lam       <- getLamb(G, type = "ETEL")
#   inner_pro <- drop(G %*% lam$lambda)
#   m         <- max(inner_pro)
#   logpi     <- inner_pro - (m + log(sum(exp(inner_pro - m))))
#   list(
#     df      = tibble(index = seq_along(y), inner_pro = inner_pro, logpi = logpi),
#     loglike = sum(logpi)
#   )
# }
 
# etel_fit_manual <- function(tet, y) {
#   objective_func <- function(theta) {
#     if (!in_hull(theta, y)) return(1e10)
#     -etel_manual(theta, y)$loglike
#   }
#   est    <- optim(tet, objective_func, method = "Nelder-Mead",
#                   control = list(maxit = 2000, reltol = 1e-10))
#   at_est <- etel_manual(est$par, y)
#   list(coef = est$par, loglike = -est$value, conv = est$convergence, df = at_est$df)
# }
 
# safe_fit <- possibly(etel_fit_manual, otherwise = NULL)
 
# biased_start <- function(yy, frac, dir, mu0 = mu, s0 = sigma2, tol = 1e-3) {
#   c_max <- max(log(yy))
#   t_hi  <- (c_max - mu0 - 1.5 * s0) / (dir[1] + 1.5 * dir[2])   # wedge edge; hull edge is inside
#   t_lo  <- 0
#   while (t_hi - t_lo > tol) {
#     t_mid <- (t_lo + t_hi) / 2
#     if (in_hull(c(mu0, s0) + t_mid * dir, yy)) t_lo <- t_mid else t_hi <- t_mid
#   }
#   list(tet0 = c(mu0, s0) + frac * t_lo * dir, t_hull = t_lo)
# }
 
 
# # -------------------- load data ----------------
 
# all_data <- read_rds("outputs/02_all_data_sims.rds")
 
 
# # -------------------- ETEL from biased start ----------------
 

 
# count(gel_results_bias, pct, delta, status, conv)
 
# write_rds(gel_results_bias, "outputs/04_gel_results_bias.rds", compress = "gz")
