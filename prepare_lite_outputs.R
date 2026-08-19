# ============================================================================
# prepare_lite_outputs.R
#
# Run this ONCE, locally, BEFORE shinylive::export(). It creates slimmed-down
# copies of your data in simulation_explorer_app/outputs/, dropping the heavy
# fit_ETEL (full "gel" S3 object) columns that your Shiny app doesn't actually
# need at runtime — only the extracted beta0_hat/beta1_hat/gt_lambda/
# implied_prob values it was already pulling out of them.
#
# This should fix the "R character strings are limited to 2^31-1 bytes" error
# from shinylive::export(), which happens when the total embedded data is
# too large for a single JSON string.
# ============================================================================

library(readr)
library(dplyr)
library(purrr)

# ---- Adjust these two paths if your folders are laid out differently -------
src_dir  <- "outputs"                              # your ORIGINAL full-size files
dest_dir <- "simulation_explorer_app/outputs"       # where the lite versions go

# Keep ALL simulation replicates.
DEPLOY_SIMS <- 1:50

# Coarsen the out_eps_mean sweep instead of dropping sims: 21 values -> 11.
# Set to NULL to keep the full seq(-5, 5, by = 0.5) sweep (no coarsening).
DEPLOY_EPS <- seq(-5, 5, by = 1)

# Thin the (th1, th2) log-likelihood grid by this factor in EACH dimension.
# 2 means "keep every 2nd grid point" -> ~4x fewer points overall.
# Set to 1 to keep the full grid resolution (no thinning).
LL_GRID_THIN <- 2

dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

apply_deploy_filters <- function(df) {
  df <- df %>% filter(sim %in% DEPLOY_SIMS)
  if (!is.null(DEPLOY_EPS)) df <- df %>% filter(out_eps_mean %in% DEPLOY_EPS)
  df
}

# ---- 1. Data + OLS results: apply sim + out_eps_mean filters -----------------
write_rds(
  read_rds(file.path(src_dir, "02_all_data_sims.rds")) %>% apply_deploy_filters(),
  file.path(dest_dir, "02_all_data_sims.rds"), compress = "gz"
)
write_rds(
  read_rds(file.path(src_dir, "02_ols_results_wide.rds")) %>% apply_deploy_filters(),
  file.path(dest_dir, "02_ols_results_wide.rds"), compress = "gz"
)

# ---- 2. GEL coefficient-only files: keep just beta0_hat/beta1_hat ----------
# (used by the Intercept/Slope density tabs)
slim_gel_coefs <- function(path_in) {
  read_rds(path_in) %>%
    apply_deploy_filters() %>%
    mutate(
      beta0_hat = map_dbl(fit_ETEL, ~ .x$coefficients[1]),
      beta1_hat = map_dbl(fit_ETEL, ~ .x$coefficients[2])
    ) %>%
    select(pct, out_eps_mean, sim, beta0_hat, beta1_hat)
}

write_rds(
  slim_gel_coefs(file.path(src_dir, "04_gel_results_bias.rds")),
  file.path(dest_dir, "05_gel_coefs_bias.rds"), compress = "gz"
)
write_rds(
  slim_gel_coefs(file.path(src_dir, "04_gel_results.rds")),
  file.path(dest_dir, "05_gel_coefs.rds"), compress = "gz"
)

# ---- 3. GEL "full" files: keep only gt_lambda + implied_prob ----------------
# (the nested x/y/is_outlier data is dropped here since it's identical to
# all_data — no need to store it twice; the app will look it up from
# all_data instead)
# (used by the gt x lambda / Implied Probabilities tabs)
slim_gel_full <- function(path_in) {
  read_rds(path_in) %>%
    apply_deploy_filters() %>%
    select(pct, out_eps_mean, sim, gt_lambda, implied_prob)
}

write_rds(
  slim_gel_full(file.path(src_dir, "04_gel_results_bias2.rds")),
  file.path(dest_dir, "05_gel_full_bias.rds"), compress = "gz"
)
write_rds(
  slim_gel_full(file.path(src_dir, "04_gel_results2.rds")),
  file.path(dest_dir, "05_gel_full.rds"), compress = "gz"
)

# ---- 4. Log-likelihood surface: filter + thin the theta grid -----------------
ll_surface_raw <- read_rds(file.path(src_dir, "06_ll_surface.rds")) %>%
  apply_deploy_filters()

if (LL_GRID_THIN > 1) {
  # Keep every Nth grid point in each dimension, based on position within
  # the sorted unique values of th1/th2 (robust to floating-point spacing).
  th1_keep <- sort(unique(ll_surface_raw$th1))[c(TRUE, rep(FALSE, LL_GRID_THIN - 1))]
  th2_keep <- sort(unique(ll_surface_raw$th2))[c(TRUE, rep(FALSE, LL_GRID_THIN - 1))]
  
  ll_surface_raw <- ll_surface_raw %>%
    filter(th1 %in% th1_keep, th2 %in% th2_keep)
}

write_rds(ll_surface_raw, file.path(dest_dir, "06_ll_surface.rds"), compress = "gz")

# ---- 5. Report file sizes so you can see what's left --------------------------
cat("\nFile sizes in", dest_dir, ":\n")
list.files(dest_dir, full.names = TRUE) %>%
  walk(~ cat(sprintf("  %-30s %6.1f MB\n", basename(.x), file.info(.x)$size / 1e6)))
