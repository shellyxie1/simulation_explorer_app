# ============================================================================
# export_sepcond.R
#
# The separation conditions of Theorem thm:compare (chapter 4 of the thesis),
# evaluated on every replicate of the linear-regression study and written to
# ../docs/linear/data/sepcond.json for the "Separation conditions" page.
# Run it from the linear/ folder after export_data.R:
#
#   Rscript export_sepcond.R
#
# Base R only, like export_data.R.
#
# For each replicate (pct, delta, sim), with n = 1000 and
#   E            = the realised outliers (is_outlier),            m = |E|
#   theta_tilde  = the displaced-start maximizer (05_gel_full_bias.rds)
#   theta        = (a) the true-start maximizer (05_gel_full.rds)
#                  (b) the truth theta_0 = (beta_0, beta_1) = (2, -3)
# the inner products eta_i(theta) = lambda_hat(theta)' g(x_i, y_i, theta) give
#   M^-_E(theta) = -max_{i in E} eta_i      M^+_E(theta) = min_{i in E} eta_i
#   b_E(theta)   =  max_{i not in E} eta_i  M'(theta_tilde) = max_i |eta_i(theta_tilde)|
# and the two conditions
#   negative:  M^-_E(theta) > (n log n + 2 n M'(theta_tilde)) / m
#   positive:  M^+_E(theta) > 0  and
#              M^+_E(theta) > b_E(theta) - log m + (n log n + 2 n M'(theta_tilde)) / (n - m).
# Besides E = the outliers, the best group of each size is tried: for the
# negative condition the k most negative inner products, for the positive
# condition the k largest, k = 1, ..., n - 1. The largest k for which the
# condition holds is exported (0 if none), with M^-_E for that group and how
# many of its members are outliers.
#
# The stored fits carry eta_i at the two maximizers. At theta_0 the inner
# problem is solved here: lambda_hat minimises mean_i exp(lambda' g_i), a
# smooth convex problem in three unknowns, by Newton's method with a
# backtracking line search on log mean exp (no overflow when eta ~ -1e5).
# ============================================================================

src  <- "outputs"
dest <- "../docs/linear/data"

THETA0 <- c(b0 = 2, b1 = -3)

# The moment conditions of the regression study (00_moments.R of the pipeline
# that produced outputs/): e, e x and e^3 with e = y - beta_0 - beta_1 x.
g_linear_reg <- function(tet, x, y) {
  e <- y - tet[1] - tet[2] * x
  cbind(e, e * x, e^3)
}

# ---- inner problem at a fixed theta -----------------------------------------
# Columns of G are rescaled for conditioning (e^3 is far larger than e for an
# outlier); eta = G lambda does not depend on the scaling.
et_eta <- function(G, tol = 1e-22, maxit = 200) {
  sc <- apply(G, 2, sd); sc[!is.finite(sc) | sc == 0] <- 1
  Gs <- sweep(G, 2, sc, "/")
  lam <- numeric(ncol(Gs))
  fval <- function(l) { e <- drop(Gs %*% l); mx <- max(e); mx + log(mean(exp(e - mx))) }
  f <- fval(lam); dec <- Inf
  for (it in seq_len(maxit)) {
    e  <- drop(Gs %*% lam); mx <- max(e); p <- exp(e - mx); p <- p / sum(p)
    gr <- drop(crossprod(Gs, p))
    H  <- crossprod(Gs * sqrt(p)) - tcrossprod(gr)
    dir <- tryCatch(-solve(H, gr), error = function(err) -solve(H + diag(1e-14, ncol(H)), gr))
    dec <- -sum(gr * dir)                       # squared Newton decrement
    if (!is.finite(dec) || dec < tol) break
    # Damped steps until the quadratic regime, then full Newton steps: the
    # extreme eta_i are sensitive to lambda, so the iteration is run to machine
    # precision rather than stopped once the objective has settled.
    step <- 1
    if (dec > 1e-6) repeat {
      fn <- fval(lam + step * dir)
      if (fn <= f - 1e-4 * step * dec || step < 1e-10) break
      step <- step / 2
    }
    lam <- lam + step * dir; f <- fval(lam)
  }
  list(eta = drop(Gs %*% lam), conv = is.finite(dec) && dec < 1e-10)
}

loglik <- function(eta) { mx <- max(eta); sum(eta - (mx + log(sum(exp(eta - mx))))) }

# ---- functionals and the best group of each size ----------------------------
functionals <- function(eta, E)
  c(Mneg = -max(eta[E]), Mpos = min(eta[E]), bE = max(eta[-E]))

groups <- function(eta, Mp, E) {
  n <- length(eta); a <- n * log(n) + 2 * n * Mp
  o <- order(eta); s <- eta[o]; k <- seq_len(n - 1)
  neg  <- -s[k] > a / k
  kneg <- if (any(neg)) max(which(neg)) else 0L
  sd   <- rev(s)                                  # descending: top k and the next one
  pos  <- sd[k] > 0 & sd[k] > sd[k + 1] - log(k) + a / (n - k)
  kpos <- if (any(pos)) max(which(pos)) else 0L
  c(kneg  = kneg,
    Mnegk = if (kneg > 0) -s[kneg] else NA_real_,
    kout  = if (kneg > 0) sum(o[seq_len(kneg)] %in% E) else NA_real_,
    kpos  = kpos)
}

# ---- one replicate -----------------------------------------------------------
one_rep <- function(x, y, E, eta_t, eta_d) {
  n <- length(y); m <- length(E)
  Mp <- max(abs(eta_d))
  f0 <- et_eta(g_linear_reg(THETA0, x, y))
  Ft <- functionals(eta_t, E); Gt <- groups(eta_t, Mp, E)
  if (f0$conv) { F0 <- functionals(f0$eta, E); G0 <- groups(f0$eta, Mp, E); ll0 <- loglik(f0$eta) }
  else { F0 <- c(Mneg = NA, Mpos = NA, bE = NA); G0 <- c(kneg = NA, Mnegk = NA, kout = NA, kpos = NA); ll0 <- NA }
  c(m = m, Mp = Mp, ll_d = loglik(eta_d), ll_t = loglik(eta_t), ll_0 = ll0, conv_0 = as.numeric(f0$conv),
    Mneg_t = Ft[["Mneg"]], Mpos_t = Ft[["Mpos"]], bE_t = Ft[["bE"]],
    kneg_t = Gt[["kneg"]], Mnegk_t = Gt[["Mnegk"]], kout_t = Gt[["kout"]], kpos_t = Gt[["kpos"]],
    Mneg_0 = F0[["Mneg"]], Mpos_0 = F0[["Mpos"]], bE_0 = F0[["bE"]],
    kneg_0 = G0[["kneg"]], Mnegk_0 = G0[["Mnegk"]], kout_0 = G0[["kout"]], kpos_0 = G0[["kpos"]])
}

# ---- run over every replicate, in the canonical slot order -------------------
message("reading ", src, "/ ...")
all_data <- readRDS(file.path(src, "02_all_data_sims.rds"))
gel      <- readRDS(file.path(src, "05_gel_full.rds"))
gel_b    <- readRDS(file.path(src, "05_gel_full_bias.rds"))

pcts <- sort(unique(all_data$pct))
epss <- sort(unique(all_data$out_eps_mean))
sims <- sort(unique(all_data$sim))

rep_key <- function(pct, eps, sim) paste(pct, eps, sim, sep = "\r")
data_by <- split(all_data, rep_key(all_data$pct, all_data$out_eps_mean, all_data$sim))
inner_by <- function(res) {
  v <- lapply(res$gt_lambda, function(g) if (is.null(g)) NULL else as.numeric(g))
  names(v) <- rep_key(res$pct, res$out_eps_mean, res$sim)
  v
}
eta_t_by <- inner_by(gel)
eta_d_by <- inner_by(gel_b)

rows <- list(); done <- 0L
for (pct in pcts) for (eps in epss) for (sim in sims) {
  k <- rep_key(pct, eps, sim)
  d <- data_by[[k]]
  et <- eta_t_by[[k]]; ed <- eta_d_by[[k]]
  rows[[k]] <- if (is.null(et) || is.null(ed) || !any(d$is_outlier) || all(d$is_outlier)) NULL
               else one_rep(d$x, d$y, which(d$is_outlier), et, ed)
  done <- done + 1L
  if (done %% 200L == 0L) message("  ", done, " / ", length(pcts) * length(epss) * length(sims))
}
fields <- c("m", "Mp", "ll_d", "ll_t", "ll_0", "conv_0",
            "Mneg_t", "Mpos_t", "bE_t", "kneg_t", "Mnegk_t", "kout_t", "kpos_t",
            "Mneg_0", "Mpos_0", "bE_0", "kneg_0", "Mnegk_0", "kout_0", "kpos_0")
tab <- t(vapply(rows, function(r) if (is.null(r)) setNames(rep(NA_real_, length(fields)), fields) else r[fields],
                numeric(length(fields))))
message(sprintf("%d replicates; inner problem at theta_0 converged in %d", nrow(tab), sum(tab[, "conv_0"] == 1, na.rm = TRUE)))

# ---- sepcond.json ------------------------------------------------------------
jnum <- function(v, digits = 7) {
  v <- as.numeric(v); s <- formatC(v, digits = digits, format = "g")
  s[!is.finite(v)] <- "null"
  paste0("[", paste(trimws(s), collapse = ","), "]")
}
jobj <- function(kv) paste0("{", paste0('"', names(kv), '":', unlist(kv), collapse = ","), "}")
out <- c(list(n = nrow(data_by[[1]]), theta0 = jnum(THETA0)), lapply(setNames(fields, fields), function(f) jnum(tab[, f])))
writeLines(jobj(out), file.path(dest, "sepcond.json"))
message("wrote ", file.path(dest, "sepcond.json"))
