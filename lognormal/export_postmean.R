# ============================================================================
# export_postmean.R
#
# The "Posterior means" tab: one BETEL posterior mean per replication for
# LN-N1000 (no outliers) and LN-N1000-C10-S3 (10% outliers, delta = 3),
# written to ../docs/lognormal/data/postmean.json. Run it from the lognormal/
# folder:
#
#   Rscript export_postmean.R
#
# Base R only, like export_data.R.
#
# The MCMC output is not in this repo, so MCMC_FILE points at it. It holds the
# mean of each chain's draws after burn-in, one row per replication -- the
# same numbers as averaging 05_MCMC_all.rds by replication, which is what the
# thesis figure LN-N1000_postmean_density.pdf does (checked equal, 2026-11-12).
# The contaminated datasets are the pct "10", delta 3 samples of the other
# tabs, so the site can link each replication to its "Simulated data" view;
# the clean datasets (pct "00") are not in the explorer.
# ============================================================================

MCMC_FILE <- path.expand(
  "~/Desktop/thesis/Check_outlier_ME_SX/check_outlier_log_normal/outputs/05_MCMC_post_means.rds"
)
dest <- "../docs/lognormal/data"

THETA_TRUE <- c(mu = 0, sigma2 = 0.25)
LABELS     <- c("00" = "LN-N1000", "10" = "LN-N1000-C10-S3")

# ---- minimal JSON writing (as in export_data.R) -----------------------------
jnum <- function(v, digits = 7) {
  s <- formatC(as.numeric(v), digits = digits, format = "g")
  s[!is.finite(as.numeric(v))] <- "null"
  paste0("[", paste(trimws(s), collapse = ","), "]")
}
jq   <- function(s) paste0('"', s, '"')
jarr <- function(items) paste0("[", paste(items, collapse = ","), "]")
jobj <- function(...) {
  kv <- list(...)
  paste0("{", paste0('"', names(kv), '":', unlist(kv), collapse = ","), "}")
}

message("reading ", MCMC_FILE, " ...")
pm <- as.data.frame(readRDS(MCMC_FILE))
pm <- pm[order(pm$pct, pm$sim), ]
stopifnot(setequal(unique(pm$pct), names(LABELS)),
          all(table(pm$pct) == 50))

scenarios <- vapply(names(LABELS), function(p) {
  d <- pm[pm$pct == p, ]
  stopifnot(length(unique(d$delta)) == 1)
  message(sprintf("  %-16s %d replications, median posterior mean mu %.3f, sigma2 %.3f",
                  LABELS[[p]], nrow(d), median(d$mu), median(d$sigma2)))
  jobj(
    label        = jq(LABELS[[p]]),
    n            = 1000,
    contaminated = tolower(p != "00"),
    pct          = jq(p),
    delta        = d$delta[1],
    sim          = jnum(d$sim),
    mu           = jnum(d$mu),
    sigma2       = jnum(d$sigma2)
  )
}, character(1))

out <- file.path(dest, "postmean.json")
writeLines(jobj(
  truth     = jobj(mu = THETA_TRUE[["mu"]], sigma2 = THETA_TRUE[["sigma2"]]),
  n         = jnum(1000),
  scenarios = jarr(scenarios)
), out)
message("wrote ", out, sprintf(" (%.1f KB)", file.size(out) / 1024))
