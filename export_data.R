# ============================================================================
# export_data.R
#
# Turns outputs/*.rds into the static payload the D3 site in docs/ reads.
# Run it from the repo root whenever the simulation outputs change:
#
#   Rscript export_data.R
#
# Base R only - no jsonlite, no tidyverse - so it runs anywhere R does.
#
# The split is deliberate. Everything the site needs *before* the reader picks
# a replicate (the factor levels, the theta grid, all 1650x3 coefficient pairs)
# goes into two small JSON files loaded once at startup. Everything that is
# per-replicate goes into one little binary per replicate, fetched on demand:
# as JSON the per-observation columns would be ~75 MB of text, and the browser
# would parse a megabyte of digits every time a slider moved.
# ============================================================================

src  <- "outputs"
dest <- "docs/data"

BETA_TRUE <- c(b0 = 2, b1 = -3)

message("reading ", src, "/ ...")
all_data     <- readRDS(file.path(src, "02_all_data_sims.rds"))
ols_wide     <- readRDS(file.path(src, "02_ols_results_wide.rds"))
gel_coefs_b  <- readRDS(file.path(src, "05_gel_coefs_bias.rds"))
gel_coefs    <- readRDS(file.path(src, "05_gel_coefs.rds"))
gel_full_b   <- readRDS(file.path(src, "05_gel_full_bias.rds"))
gel_full     <- readRDS(file.path(src, "05_gel_full.rds"))
ll_surface   <- readRDS(file.path(src, "06_ll_surface.rds"))

# ---- factor levels, in the order every array below is written --------------
pcts <- sort(unique(all_data$pct))
epss <- sort(unique(all_data$out_eps_mean))
sims <- sort(unique(all_data$sim))
th1  <- sort(unique(ll_surface$th1))
th2  <- sort(unique(ll_surface$th2))

n_obs <- nrow(all_data) / (length(pcts) * length(epss) * length(sims))
stopifnot(n_obs == round(n_obs))
n_obs <- as.integer(n_obs)

message(sprintf(
  "  %d pct x %d eps x %d sim = %d replicates of %d observations",
  length(pcts), length(epss), length(sims),
  length(pcts) * length(epss) * length(sims), n_obs
))

# ---- minimal JSON writing ---------------------------------------------------
# Only what this script emits: numeric vectors, character vectors and nested
# objects built from them. `format = "g"` keeps 2 as "2" rather than "2.000000"
# and drops to exponent notation only where it actually saves characters.
jnum <- function(v, digits = 7) {
  s <- formatC(as.numeric(v), digits = digits, format = "g")
  s[!is.finite(as.numeric(v))] <- "null"
  paste0("[", paste(trimws(s), collapse = ","), "]")
}
jstr <- function(v) paste0("[", paste0('"', v, '"', collapse = ","), "]")
jobj <- function(...) {
  kv <- list(...)
  paste0("{", paste0('"', names(kv), '":', unlist(kv), collapse = ","), "}")
}

dir.create(file.path(dest, "rep"), recursive = TRUE, showWarnings = FALSE)

# ---- meta.json --------------------------------------------------------------
# `layout` mirrors the byte layout written below. The site reads its offsets
# from here rather than hard-coding them, so changing n_obs or thinning the
# theta grid needs no edit on the JavaScript side.
n_ll <- length(th1) * length(th2)
writeLines(jobj(
  n_obs     = n_obs,
  beta_true = jobj(b0 = BETA_TRUE[["b0"]], b1 = BETA_TRUE[["b1"]]),
  pct       = jstr(pcts),
  eps       = jnum(epss),
  # The filename fragment for each eps, so the site never has to reproduce R's
  # number formatting to build a URL.
  eps_slug  = jstr(formatC(epss, format = "g")),
  sim       = jnum(sims),
  th1       = jnum(th1),
  th2       = jnum(th2),
  # Byte offsets into a replicate file. `ll` is row-major over th1, so
  # ll[i1 * length(th2) + i2] is the log-likelihood at (th1[i1], th2[i2]).
  layout    = jobj(
    x       = 0L,
    y       = 4L * n_obs,
    gt_bias = 8L * n_obs,
    gt      = 12L * n_obs,
    ll      = 16L * n_obs,
    outlier = 16L * n_obs + 4L * n_ll,
    bytes   = 17L * n_obs + 4L * n_ll
  )
), file.path(dest, "meta.json"))

# ---- coefs.json -------------------------------------------------------------
# One (beta0_hat, beta1_hat) pair per method per replicate, in the canonical
# order pct-major then eps then sim, so the site indexes it arithmetically
# instead of carrying 1650 copies of the three key columns.
key <- expand.grid(sim = sims, out_eps_mean = epss, pct = pcts,
                   stringsAsFactors = FALSE)
key$slot <- seq_len(nrow(key))

align <- function(df) {
  m <- merge(key, df, by = c("pct", "out_eps_mean", "sim"), all.x = TRUE,
             sort = FALSE)
  m <- m[order(m$slot), ]
  stopifnot(nrow(m) == nrow(key), !anyNA(m$beta0_hat))
  m
}

coef_block <- function(df) {
  m <- align(df)
  jobj(b0 = jnum(m$beta0_hat), b1 = jnum(m$beta1_hat))
}

writeLines(jobj(
  OLS      = coef_block(ols_wide),
  gel_bias = coef_block(gel_coefs_b),
  gel      = coef_block(gel_coefs)
), file.path(dest, "coefs.json"))

# ---- per-replicate binaries -------------------------------------------------
# Layout, little-endian, all float32 except the trailing flags:
#
#   x[n_obs]  y[n_obs]  gt_bias[n_obs]  gt[n_obs]  ll[n_ll]  outlier[n_obs]:u8
#
# Only gt = gt %*% lambda is stored for each GEL fit: the implied probabilities
# are exactly softmax(gt) (checked to 1e-15 against the fitted values), so the
# site derives them and the payload halves.
#
# The three tables are split once by replicate key rather than filtered 1650
# times each; the repeated filter was the whole runtime of an earlier draft.
rep_key <- function(pct, eps, sim) paste(pct, eps, sim, sep = "\r")
slug    <- function(pct, eps, sim) sprintf("%s-%s-%s", pct, formatC(eps, format = "g"), sim)

split_by_key <- function(df) split(df, rep_key(df$pct, df$out_eps_mean, df$sim))

message("splitting tables ...")
data_by  <- split_by_key(all_data)
fullb_by <- split_by_key(gel_full_b)
full_by  <- split_by_key(gel_full)
ll_by    <- split_by_key(ll_surface)

message("writing ", length(data_by), " replicate files ...")
written <- 0L
for (pct in pcts) for (eps in epss) for (sim in sims) {
  k <- rep_key(pct, eps, sim)
  d <- data_by[[k]]
  ll <- ll_by[[k]]
  stopifnot(nrow(d) == n_obs, nrow(ll) == n_ll)

  # The theta grid is written row-major over th1; `match` rather than an
  # assumed row order, since the rds is only grouped, not sorted.
  z <- numeric(n_ll)
  z[(match(ll$th1, th1) - 1L) * length(th2) + match(ll$th2, th2)] <- ll$ll
  stopifnot(!anyNA(z))

  con <- file(file.path(dest, "rep", paste0(slug(pct, eps, sim), ".bin")), "wb")
  writeBin(as.numeric(d$x),                          con, size = 4, endian = "little")
  writeBin(as.numeric(d$y),                          con, size = 4, endian = "little")
  writeBin(as.numeric(fullb_by[[k]]$gt_lambda[[1]]), con, size = 4, endian = "little")
  writeBin(as.numeric(full_by[[k]]$gt_lambda[[1]]),  con, size = 4, endian = "little")
  writeBin(z,                                        con, size = 4, endian = "little")
  writeBin(as.integer(d$is_outlier),                 con, size = 1)
  close(con)

  written <- written + 1L
  if (written %% 200L == 0L) message("  ", written, " / ", length(data_by))
}

sz <- sum(file.info(list.files(dest, recursive = TRUE, full.names = TRUE))$size)
message(sprintf("done - %s holds %.1f MB across %d files",
                dest, sz / 1e6,
                length(list.files(dest, recursive = TRUE))))
