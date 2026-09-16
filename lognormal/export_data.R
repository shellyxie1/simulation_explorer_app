# ============================================================================
# export_data.R
#
# Turns outputs/*.rds into the static payload the D3 site in docs/ reads.
# Run it from the lognormal/ folder whenever the simulation outputs change:
#
#   Rscript export_data.R
#
# Base R only - no jsonlite, no tidyverse - so it runs anywhere R does.
#
# The split follows the linear-regression explorer this was adapted from.
# Everything the site needs *before* the reader picks a replicate (the factor
# levels, the (mu, sigma2) grid, all 1200 x 2 estimates) goes into two small
# JSON files loaded once at startup. Everything that is per-replicate goes into
# one little binary per replicate, fetched on demand.
#
# Differences from the linear version:
#   * design is (pct, delta, sim); parameters are (mu, sigma2); no OLS
#   * there is no x column: y is indexed by observation number
#   * the log-likelihood grid has infeasible cells (outside the convex hull of
#     the moment conditions). They go out as NaN in the float32 block plus an
#     explicit hull flag, so the site can shade them.
#   * the biased start point of each ETEL_bias fit is exported too, so the
#     surface view can show where the optimiser began.
# ============================================================================

src  <- "outputs"
dest <- "../docs/lognormal/data"

THETA_TRUE <- c(mu = 0, sigma2 = 0.25)

message("reading ", src, "/ ...")
all_data   <- readRDS(file.path(src, "02_all_data_sims.rds"))
gel_b      <- readRDS(file.path(src, "04_gel_results_bias.rds"))
gel        <- readRDS(file.path(src, "04_gel_results.rds"))
ll_surface <- readRDS(file.path(src, "06_ll_surface.rds"))

# ---- factor levels, in the order every array below is written --------------
pcts   <- sort(unique(all_data$pct))        # character: "01" "05" "10"
deltas <- sort(unique(all_data$delta))      # numeric:   -4 ... 4 (no 0)
sims   <- sort(unique(all_data$sim))
mus    <- sort(unique(ll_surface$mu))
s2s    <- sort(unique(ll_surface$sigma2))

n_rep <- length(pcts) * length(deltas) * length(sims)
n_obs <- nrow(all_data) / n_rep
stopifnot(n_obs == round(n_obs))
n_obs <- as.integer(n_obs)
n_ll  <- length(mus) * length(s2s)

message(sprintf(
  "  %d pct x %d delta x %d sim = %d replicates of %d observations; %d x %d ll grid",
  length(pcts), length(deltas), length(sims), n_rep, n_obs, length(mus), length(s2s)
))

# ---- minimal JSON writing ---------------------------------------------------
# Only what this script emits: numeric vectors, character vectors and nested
# objects built from them. NA / NaN / Inf become null, which is the only way
# JSON can carry a missing number.
jnum <- function(v, digits = 7) {
  v <- as.numeric(v)
  s <- formatC(v, digits = digits, format = "g")
  s[!is.finite(v)] <- "null"
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
# from here rather than hard-coding them, so changing n_obs or the theta grid
# needs no edit on the JavaScript side.
writeLines(jobj(
  n_obs      = n_obs,
  theta_true = jobj(mu = THETA_TRUE[["mu"]], sigma2 = THETA_TRUE[["sigma2"]]),
  pct        = jstr(pcts),
  delta      = jnum(deltas),
  # The filename fragment for each delta, so the site never has to reproduce
  # R's number formatting to build a URL.
  delta_slug = jstr(formatC(deltas, format = "g")),
  sim        = jnum(sims),
  mu         = jnum(mus),
  sigma2     = jnum(s2s),
  # Byte offsets into a replicate file. `ll` and `hull` are row-major over mu,
  # so ll[i * length(sigma2) + j] is the log-likelihood at (mu[i], sigma2[j]).
  layout     = jobj(
    y       = 0L,
    gt_bias = 4L * n_obs,
    gt      = 8L * n_obs,
    ll      = 12L * n_obs,
    hull    = 12L * n_obs + 4L * n_ll,
    outlier = 12L * n_obs + 5L * n_ll,
    bytes   = 13L * n_obs + 5L * n_ll
  )
), file.path(dest, "meta.json"))

# ---- coefs.json -------------------------------------------------------------
# One (mu_hat, sigma2_hat, loglike, conv) per method per replicate, in the
# canonical order pct-major then delta then sim, so the site indexes it
# arithmetically instead of carrying 1200 copies of the three key columns.
key <- expand.grid(sim = sims, delta = deltas, pct = pcts, stringsAsFactors = FALSE)
key$slot <- seq_len(nrow(key))

pick <- function(lst, k) vapply(lst, function(v) if (is.null(v)) NA_real_ else v[[k]], 0)

scalar_table <- function(res, biased = FALSE) {
  df <- data.frame(
    pct        = res$pct,
    delta      = res$delta,
    sim        = res$sim,
    mu_hat     = pick(res$coef, 1),
    sigma2_hat = pick(res$coef, 2),
    loglike    = as.numeric(res$loglike),
    conv       = as.numeric(res$conv),
    stringsAsFactors = FALSE
  )
  if (biased) {
    df$start_mu     <- pick(res$tet0, 1)
    df$start_sigma2 <- pick(res$tet0, 2)
    df$t_hull       <- as.numeric(res$t_hull)
  }
  df
}

align <- function(df) {
  m <- merge(key, df, by = c("pct", "delta", "sim"), all.x = TRUE, sort = FALSE)
  m <- m[order(m$slot), ]
  stopifnot(nrow(m) == nrow(key))
  m
}

coef_block <- function(df, biased = FALSE) {
  m <- align(df)
  if (biased) {
    jobj(mu = jnum(m$mu_hat), sigma2 = jnum(m$sigma2_hat),
         ll = jnum(m$loglike), conv = jnum(m$conv),
         start_mu = jnum(m$start_mu), start_sigma2 = jnum(m$start_sigma2),
         t_hull = jnum(m$t_hull))
  } else {
    jobj(mu = jnum(m$mu_hat), sigma2 = jnum(m$sigma2_hat),
         ll = jnum(m$loglike), conv = jnum(m$conv))
  }
}

writeLines(jobj(
  gel_bias = coef_block(scalar_table(gel_b, biased = TRUE), biased = TRUE),
  gel      = coef_block(scalar_table(gel))
), file.path(dest, "coefs.json"))

# ---- per-replicate binaries -------------------------------------------------
# Layout, little-endian, float32 blocks first and the byte flags last:
#
#   y[n_obs]  gt_bias[n_obs]  gt[n_obs]  ll[n_ll]  hull[n_ll]:u8  outlier[n_obs]:u8
#
# Only the inner product lambda' g_i is stored for each ETEL fit: the implied
# probabilities are exactly softmax of it (log pi_i = inner_pro_i - logsumexp,
# which is how 00_helper_func.R computes them), so the site derives them.
# Infeasible grid cells are NaN in `ll` and 0 in `hull`.
rep_key <- function(pct, delta, sim) paste(pct, delta, sim, sep = "\r")
slug    <- function(pct, delta, sim)
  sprintf("%s-%s-%s", pct, formatC(delta, format = "g"), sim)

message("splitting tables ...")
data_by <- split(all_data, rep_key(all_data$pct, all_data$delta, all_data$sim))
ll_by   <- split(ll_surface, rep_key(ll_surface$pct, ll_surface$delta, ll_surface$sim))

inner_by <- function(res) {
  v <- lapply(res$fit_ETEL, function(f)
    if (is.null(f)) rep(NaN, n_obs) else as.numeric(f$df$inner_pro))
  names(v) <- rep_key(res$pct, res$delta, res$sim)
  v
}
gtb_by <- inner_by(gel_b)
gt_by  <- inner_by(gel)

message("writing ", n_rep, " replicate files ...")
written <- 0L
for (pct in pcts) for (delta in deltas) for (sim in sims) {
  k  <- rep_key(pct, delta, sim)
  d  <- data_by[[k]]
  ll <- ll_by[[k]]
  stopifnot(nrow(d) == n_obs, nrow(ll) == n_ll)
  d <- d[order(d$index), ]

  # Row-major over mu; `match` rather than an assumed row order.
  pos <- (match(ll$mu, mus) - 1L) * length(s2s) + match(ll$sigma2, s2s)
  z <- rep(NaN, n_ll); z[pos] <- ifelse(is.na(ll$ll), NaN, ll$ll)
  h <- integer(n_ll);  h[pos] <- as.integer(ll$hull)

  con <- file(file.path(dest, "rep", paste0(slug(pct, delta, sim), ".bin")), "wb")
  writeBin(as.numeric(d$y),        con, size = 4, endian = "little")
  writeBin(gtb_by[[k]],            con, size = 4, endian = "little")
  writeBin(gt_by[[k]],             con, size = 4, endian = "little")
  writeBin(z,                      con, size = 4, endian = "little")
  writeBin(h,                      con, size = 1)
  writeBin(as.integer(d$is_out),   con, size = 1)
  close(con)

  written <- written + 1L
  if (written %% 200L == 0L) message("  ", written, " / ", n_rep)
}

sz <- sum(file.info(list.files(dest, recursive = TRUE, full.names = TRUE))$size)
message(sprintf("done - %s holds %.1f MB across %d files",
                dest, sz / 1e6, length(list.files(dest, recursive = TRUE))))
