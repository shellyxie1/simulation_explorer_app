# helpers 
library(readr)
library(gmm)
library(tidyverse)
# for plotting lognormal densities
app_theme <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text   = element_text(family = "serif", colour = "black"),
      axis.title  = element_text(colour = "black"),
      axis.ticks  = element_line(colour = "black"),
      plot.title  = element_text(face = "bold"),
      legend.position = "bottom"
    )
}




# ------------------- for 03_ETEL_bias.R and 03_ETEL.R
# function to check whether convex hull satisfies
in_hull <- function(tet, y) {
  G   <- g_lognormal(tet, y)
  lp  <- lpSolve::lp(
    direction = "min", objective.in = rep(0, nrow(G)),
    const.mat = rbind(t(G), rep(1, nrow(G))),
    const.dir = rep("=", ncol(G) + 1),
    const.rhs = c(rep(0, ncol(G)), 1)
  )
  lp$status == 0
}

etel_manual <- function(tet, y) {
  G         <- g_lognormal(tet, y)
  lam       <- getLamb(G, type = "ETEL")
  inner_pro <- drop(G %*% lam$lambda)
  m         <- max(inner_pro)
  logpi     <- inner_pro - (m + log(sum(exp(inner_pro - m))))
  list(
    df      = tibble(index = seq_along(y), inner_pro = inner_pro, logpi = logpi),
    loglike = sum(logpi)
  )
}

# to avoid underflow
etel_fit_manual <- function(tet, y) {
  objective_func <- function(theta) {
    if (!in_hull(theta, y)) return(1e10)
    -etel_manual(theta, y)$loglike
  }
  est    <- optim(tet, objective_func, method = "Nelder-Mead",
                  control = list(maxit = 2000, reltol = 1e-10))
  at_est <- etel_manual(est$par, y)
  list(coef = est$par, loglike = -est$value, conv = est$convergence, df = at_est$df)
}

#if etel_fit_manual throws an error, returns NULL
safe_fit <- possibly(etel_fit_manual, otherwise = NULL)

# find biased starting values that work
# I want smt still in the convex hull

biased_start <- function(yy, frac = 0.7, dir = c(1, 1), mu0 = mu, s0 = sigma2, tol = 1e-3) {
  c_max  <- max(log(yy))
  t_hi   <- (c_max - mu0 - 1.5 * s0) / (dir[1] + 1.5 * dir[2])  # wedge edge: hull edge is inside this
  t_lo   <- 0
  while (t_hi - t_lo > tol) {
    t_mid <- (t_lo + t_hi) / 2
    if (in_hull(c(mu0, s0) + t_mid * dir, yy)) t_lo <- t_mid 
    else t_hi <- t_mid
  }
  list(tet0 = c(mu0, s0) + frac * t_lo * dir, t_hull = t_lo)
}



# ------ for 06_ll_surface.R
compute_ll_surface <- function(y) {
  theta_grid |>
    rowwise() |>
    mutate(
      hull = in_hull(c(mu, sigma2), y),
      ll   = if (!hull) NA_real_ else tryCatch(
        etel_manual(c(mu, sigma2), y)$loglike,
        error = \(e) NA_real_
      )
    ) |>
    ungroup()
}






