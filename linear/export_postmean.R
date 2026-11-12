# ============================================================================
# export_postmean.R
#
# The "Posterior means" tab: one BETEL posterior mean per replication for the
# Chapter 3 scenarios, written to ../docs/linear/data/postmean.json.
# Run it from the linear/ folder:
#
#   Rscript export_postmean.R
#
# Base R only, like export_data.R.
#
# The MCMC draws are not in this repo (twelve files of 12 MB each), so
# MCMC_DIR points at them. Each MCMC_b{0,1}_<design>.rds is a 30000 x 50
# matrix of draws after burn-in, one column per replication, and the posterior
# mean is the column mean -- the same numbers as the thesis figure
# LR-N1000_postmean_density.pdf. Only the 600 posterior means travel to the
# site; the kernel densities are drawn in the browser.
#
# These are the Chapter 3 scenarios (outlier errors from N(5, 1)), a separate
# simulation from the sweep in outputs/ that the other tabs read.
# ============================================================================

MCMC_DIR <- path.expand("~/Desktop/thesis/outputs/linear_reg/03_MCMC")
dest     <- "../docs/linear/data"

DESIGNS <- data.frame(
  design       = c("D11", "D12", "D21", "D22", "D31", "D32"),
  label        = c("LR-N100", "LR-N100-C10", "LR-N500", "LR-N500-C10",
                   "LR-N1000", "LR-N1000-C10"),
  n            = c(100, 100, 500, 500, 1000, 1000),
  contaminated = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
)
BETA_TRUE <- c(b0 = 2, b1 = -3)

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

post_mean <- function(par, design) {
  draws <- readRDS(file.path(MCMC_DIR, sprintf("MCMC_%s_%s.rds", par, design)))
  stopifnot(is.matrix(draws), ncol(draws) == 50)
  colMeans(draws)
}

message("reading ", MCMC_DIR, "/ ...")
scenarios <- vapply(seq_len(nrow(DESIGNS)), function(i) {
  d  <- DESIGNS[i, ]
  b0 <- post_mean("b0", d$design)
  b1 <- post_mean("b1", d$design)
  message(sprintf("  %-13s %d replications, median posterior mean b0 %.3f, b1 %.3f",
                  d$label, length(b0), median(b0), median(b1)))
  jobj(
    label        = jq(d$label),
    design       = jq(d$design),
    n            = d$n,
    contaminated = tolower(d$contaminated),
    sim          = jnum(seq_along(b0)),
    b0           = jnum(b0),
    b1           = jnum(b1)
  )
}, character(1))

out <- file.path(dest, "postmean.json")
writeLines(jobj(
  truth     = jobj(b0 = BETA_TRUE[["b0"]], b1 = BETA_TRUE[["b1"]]),
  n         = jnum(sort(unique(DESIGNS$n))),
  scenarios = jarr(scenarios)
), out)
message("wrote ", out, sprintf(" (%.1f KB)", file.size(out) / 1024))
