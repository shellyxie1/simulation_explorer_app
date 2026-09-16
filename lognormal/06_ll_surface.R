library(tidyverse)
library(ggplot2)
library(gmm)


source("00_moments.R")
source("00_helper_func.R")   


theta_grid <- crossing(
  mu     = seq(-1.5, 3.5, by = 0.1),
  sigma2 = seq(0.05, 2.0, by = 0.05)
)

# test <- theta_grid |>
#   rowwise()
# all_data <- read_rds("outputs/02_all_data_sims.rds")
# test_df <- all_data |>
#   filter(pct == "10", delta == 4, sim == 5)

# hull <- in_hull(test[10,],test_df$y)


# --- Main loop ----------------------------------------------------------------

all_data <- read_rds("outputs/02_all_data_sims.rds")

# runs <- crossing(
#   pct_val   = c("01", "05", "10"),
#   delta_val = c(-4, -3, -2, -1, 1, 2, 3, 4),
#   sim_val   = 1:50
# )

# all_results <- runs |>
#   pmap(function(pct_val, delta_val, sim_val) {

#     y <- all_data |>
#       filter(pct == pct_val, delta == delta_val, sim == sim_val) |>
#       pull(y)

#     compute_ll_surface(y) |>
#       mutate(pct = pct_val, delta = delta_val, sim = sim_val)
#   }) |>
#   list_rbind()

n_datasets <- nrow(distinct(all_data, pct, delta, sim))
pb <- cli::cli_progress_bar(
  "ETEL surfaces", total = n_datasets,
  format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | elapsed {cli::pb_elapsed} | ETA {cli::pb_eta}"
)

all_results <- all_data |>
  group_by(pct, delta, sim) |> # use the data
  group_modify(\(d, key) { # then use the theta_grid to compute loglikelihood
    cli::cli_progress_update(id = pb)
    compute_ll_surface(d$y)
  }) |>
  ungroup()

cli::cli_progress_done(id = pb)

# test_if_same <- all_results |>
#   relocate(pct, delta, sim, .after = last_col())

# identical(test_if_same, read_rds("outputs/06_ll_surface.rds"))
write_rds(all_results, "outputs/06_ll_surface.rds", compress = "gz")

