# plot (1- epsilon) LN(mu, sigma2) + epsilon LN(mu + delta, sigma2)

library(tidyverse)
library(gmm)
library(ggplot2)
library(readr)

# all_data <- read_rds("outputs/02_all_data_sims.rds")

source("00_helper_func.R")
# --------------------------------------------------
# s0 <- 0.25
# delta <- 3
# mu0 <- 0


plot_lognormal_dens <- function(mu0 = 0, s0 = 0.25, delta = 3){
  sd0 <- sqrt(s0)
  x <- seq(0, exp(max(mu0, mu0 + delta) + 3 * sd0), length.out = 500)
  dens <- tibble(x = x, Clean = dlnorm(x, mu0, sd0), Outlier = dlnorm(x,mu0 + delta, sd0))
  long_des <- dens |>
    pivot_longer(-x, names_to = "component", values_to = "density")
  ggplot(long_des, aes(x = x, y = density, colour = component)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values= c(Clean = "steelblue4", Outlier = "red")) +
    labs(x = "y", y = "Density", color = NULL) +
    app_theme()
}

# plot_lognormal_dens(delta = 1)










# plot_components <- function(mu0 = 0, s0 = 0.25, delta = 3, eps = NULL, log_scale = FALSE) {
#   sd0 <- sqrt(s0)
# x <- if (log_scale) seq(min(mu0, mu0 + delta) - 4 * sd0, max(mu0, mu0 + delta) + 4 * sd0, length.out = 500)
#      else            seq(0, exp(max(mu0, mu0 + delta) + 3 * sd0), length.out = 500)
#   d   <- if (log_scale) \(x, m) dnorm(x, m, sd0) else \(x, m) dlnorm(x, m, sd0)

#   dens <- tibble(x = x, Clean = d(x, mu0), Outlier = d(x, mu0 + delta))
#   if (!is.null(eps)) dens <- dens |> mutate(Mixture = (1 - eps) * Normal + eps * Outlier)

#   dens |>
#     pivot_longer(-x, names_to = "component", values_to = "density") |>
#     ggplot(aes(x = x, y = density, colour = component)) +
#     geom_line(linewidth = 1) +
#     scale_colour_manual(values = c(Clean = "steelblue4", Outlier = "red3", Mixture = "black")) +
#     labs(x = if (log_scale) "log y" else "y", y = "Density", colour = NULL,
#          title = TeX(sprintf("LN(%s, %s) vs LN(%s + %s, %s)", mu0, s0, mu0, delta, s0))) +
#     theme_bw(base_size = 14)
# }

# plot_components(delta = 4)                          # two lognormals, raw scale
# plot_components(delta = 3, log_scale = TRUE)        # two normals in log y
# plot_components(delta = 3, eps = 0.10)              # plus the mixture the data are drawn from

# plot_components(delta = -3)

# plot_components <- function(mu0 = 0, s0 = 0.25, delta = 3, eps = NULL, log_scale = FALSE,
#                             xlim = NULL) {
#   sd0 <- sqrt(s0)
#   if (is.null(xlim)) {
#     xlim <- if (log_scale) c(min(mu0, mu0 + delta) - 4 * sd0, max(mu0, mu0 + delta) + 4 * sd0)
#             else            c(0, exp(max(mu0, mu0 + delta) + 3 * sd0))
#   }
#   x <- seq(xlim[1], xlim[2], length.out = 500)
#   d <- if (log_scale) \(x, m) dnorm(x, m, sd0) else \(x, m) dlnorm(x, m, sd0)

#   dens <- tibble(x = x, Clean = d(x, mu0), Outlier = d(x, mu0 + delta))
#   if (!is.null(eps)) dens <- dens |> mutate(Mixture = (1 - eps) * Clean + eps * Outlier)

#   dens |>
#     pivot_longer(-x, names_to = "component", values_to = "density") |>
#     ggplot(aes(x = x, y = density, colour = component)) +
#     geom_line(linewidth = 1) +
#     scale_colour_manual(values = c(Clean = "steelblue4", Outlier = "red3", Mixture = "black")) +
#     coord_cartesian(xlim = xlim) +
#     labs(x = if (log_scale) "log y" else "y", y = "Density", colour = NULL,
#          title = TeX(sprintf("LN(%s, %s) vs LN(%s %+d, %s)", mu0, s0, mu0, delta, s0))) +
#     theme_bw(base_size = 14)
# }

# library(patchwork)
# xr <- c(0, 60)
# plot_components(delta =  3, xlim = xr) / plot_components(delta = -3, xlim = xr)
