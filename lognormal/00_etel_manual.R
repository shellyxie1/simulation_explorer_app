in_hull <- function(tet, y) {
  G  <- g_lognormal(tet, y)
  lp <- lpSolve::lp(
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
 
safe_fit <- possibly(etel_fit_manual, otherwise = NULL)
 
biased_start <- function(yy, frac, dir, mu0, s0, tol = 1e-3) {
  c_max <- max(log(yy))
  t_hi  <- (c_max - mu0 - 1.5 * s0) / (dir[1] + 1.5 * dir[2])
  t_lo  <- 0
  while (t_hi - t_lo > tol) {
    t_mid <- (t_lo + t_hi) / 2
    if (in_hull(c(mu0, s0) + t_mid * dir, yy)) t_lo <- t_mid else t_hi <- t_mid
  }
  list(tet0 = c(mu0, s0) + frac * t_lo * dir, t_hull = t_lo)
}