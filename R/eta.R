#' Quadratic metric used to score a fused estimate
#'
#' Returns `A = -M`, the (negative) trial-information metric of the Taylor
#' surrogate. With `target = NULL`, `M` is the full trial information, which
#' scores all parameters. With a `target` set, `M` is the *marginal* trial
#' information of those parameters (the inverse of the corresponding block of
#' the trial covariance), embedded in a zero matrix, so borrowing is judged
#' only by its effect on the targeted parameters and nuisance parameters do
#' not count against it.
#'
#' @param fit_exp the trial model.
#' @param target `NULL`, parameter names, or integer positions.
#' @param robust use the inverse sandwich covariance as the information.
#' @return a `p x p` matrix.
#' @export
mdpl_metric <- function(fit_exp, target = NULL, robust = FALSE) {
  OI <- model_info(fit_exp, robust)
  if (is.null(target)) return(-OI)
  ix <- if (is.character(target)) match(target, names(fit_exp$par)) else as.integer(target)
  if (anyNA(ix) || any(ix < 1L | ix > nrow(OI))) stop("`target` does not match the parameters")
  M <- matrix(0, nrow(OI), ncol(OI), dimnames = dimnames(OI))
  B <- solve(solve(OI)[ix, ix, drop = FALSE])
  M[ix, ix] <- (B + t(B)) / 2
  -M
}

#' Fast Taylor-surrogate selection of the borrowing weights
#'
#' Chooses the power-likelihood weights \eqn{\eta = (\eta_1, \ldots, \eta_K)
#' \in [0,1]^K}, one per observational source, by maximising the local
#' (Taylor) surrogate of the trial's predictive score,
#' \deqn{f_\lambda(\eta) = g(\eta)^\top S(\eta)^{-1} A\, S(\eta)^{-1} g(\eta) +
#'   \lambda\, \mathrm{tr}\{A\, S(\eta)^{-1}\},}
#' where \eqn{S(\eta) = I_e + \sum_k \eta_k I_k},
#' \eqn{g(\eta) = \sum_k \eta_k I_k(\hat\theta_k - \hat\theta_e)}, and
#' \eqn{A = -M} is the metric from [mdpl_metric()]. The trace multiplier
#' \eqn{\lambda} is `2` for the leave-one-out ELPD criterion and `1` for the
#' more conservative Jensen lower bound. The exact gradient is passed to
#' L-BFGS-B, and one call returns the weights for every trace multiplier.
#'
#' @param fit_exp the trial model.
#' @param fit_obs one model or a list of them (one weight per source).
#' @param target `NULL` (score all parameters) or the parameters that matter,
#'   e.g. `treatment_terms(fit_exp, "A")`. Ignored if `A` is given.
#' @param A optional `p x p` metric overriding `target`.
#' @param trace_mult numeric vector of trace multipliers, default `c(1, 2)`.
#' @param gradient use the analytic gradient (default); `FALSE` uses finite
#'   differences, for checking.
#' @param robust use the inverse sandwich covariance as the information.
#' @return an object of class `"taylor_eta"` with `eta` (a `K x
#'   length(trace_mult)` matrix, columns `t1`, `t2`, ...), `value`, `term1`,
#'   `term2`, `convergence`, `trace_mult`, `A`.
#' @export
taylor_eta <- function(fit_exp, fit_obs, target = NULL, A = NULL,
                       trace_mult = c(1, 2), gradient = TRUE, robust = FALSE) {
  if (inherits(fit_obs, "mdpl_model")) fit_obs <- list(fit_obs)
  K <- length(fit_obs)
  if (K < 1L) stop("At least one observational fit is required")
  p <- length(fit_exp$par)
  for (f in fit_obs)
    if (!identical(names(f$par), names(fit_exp$par)))
      stop("Trial and observational models must have the same parameters")
  if (is.null(A)) A <- mdpl_metric(fit_exp, target, robust)
  A <- unname(A)
  Ie <- unname(model_info(fit_exp, robust))
  phi_e <- unname(fit_exp$par)

  Is <- lapply(fit_obs, function(f) unname(model_info(f, robust)))
  Ivec <- vapply(Is, as.numeric, numeric(p * p))
  D <- vapply(seq_len(K), function(k) as.numeric(Is[[k]] %*% (fit_obs[[k]]$par - phi_e)), numeric(p))
  dim(Ivec) <- c(p * p, K)
  dim(D) <- c(p, K)

  inv_pd <- function(S) {
    ch <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(ch)) solve(S) else chol2inv(ch)
  }
  parts <- function(x) {
    Sinv <- inv_pd(Ie + matrix(Ivec %*% x, p, p))
    s <- Sinv %*% (D %*% x)
    As <- A %*% s
    list(Sinv = Sinv, s = s, As = As, term1 = sum(s * As), term2 = sum(A * Sinv))
  }
  gr_parts <- function(pt) {
    w <- pt$Sinv %*% pt$As
    B <- pt$Sinv %*% A %*% pt$Sinv
    list(g1 = 2 * (crossprod(D, w) - crossprod(Ivec, as.numeric(tcrossprod(w, pt$s)))),
         g2 = -crossprod(Ivec, as.numeric(B)))
  }

  nT <- length(trace_mult)
  eta <- matrix(NA_real_, K, nT,
                dimnames = list(names(fit_obs) %||% paste0("source", seq_len(K)),
                                paste0("t", trace_mult)))
  out <- list(value = numeric(nT), term1 = numeric(nT), term2 = numeric(nT),
              convergence = integer(nT), counts = matrix(NA_integer_, 2L, nT))
  for (j in seq_len(nT)) {
    lam <- trace_mult[j]
    fn <- function(x) { pt <- parts(x); pt$term1 + lam * pt$term2 }
    gr <- if (gradient) function(x) {
      g <- gr_parts(parts(x))
      as.numeric(g$g1 + lam * g$g2)
    }
    res <- stats::optim(par = rep(0, K), fn = fn, gr = gr, method = "L-BFGS-B",
                        lower = rep(0, K), upper = rep(1, K),
                        control = list(maxit = 1000, factr = 1e-7, pgtol = 1e-8,
                                       fnscale = -1))
    pt <- parts(res$par)
    eta[, j] <- res$par
    out$value[j] <- res$value
    out$term1[j] <- pt$term1
    out$term2[j] <- pt$term2
    out$convergence[j] <- res$convergence
    out$counts[, j] <- res$counts
  }
  structure(c(list(eta = eta), out, list(trace_mult = trace_mult, A = A, K = K)),
            class = "taylor_eta")
}

#' @export
print.taylor_eta <- function(x, digits = 3, ...) {
  cat("Taylor-surrogate borrowing weights (K =", x$K, "sources)\n")
  print(round(x$eta, digits))
  if (any(x$convergence != 0)) warning("L-BFGS-B did not converge for some trace multipliers")
  invisible(x)
}

#' Select the borrowing weights by the trial's expected log predictive density
#'
#' For each candidate weight, draws from the normal approximation of the power
#' posterior and scores the **trial** data by WAIC (log pointwise predictive
#' density minus the effective number of parameters). With one source the
#' weight is chosen on a grid; with several sources the grid is searched one
#' coordinate at a time (two sweeps). The same standard normal draws are used
#' for every candidate, which makes the comparison between weights much less
#' noisy.
#'
#' @param fit_exp trial model.
#' @param fit_obs one model or a list of them.
#' @param grid candidate weights for each source.
#' @param n_draws posterior draws per candidate.
#' @param robust use the inverse sandwich covariance as the information.
#' @return an object of class `"elpd_eta"` with `eta` (a `K x 1` matrix,
#'   column `elpd`) and `table` (candidate weight and ELPD for the last
#'   coordinate sweep, source by source in `tables`).
#' @export
elpd_eta <- function(fit_exp, fit_obs, grid = c(seq(0, 0.2, 0.04), seq(0.25, 1, 0.05)),
                     n_draws = 1000L, robust = FALSE) {
  if (inherits(fit_obs, "mdpl_model")) fit_obs <- list(fit_obs)
  if (is.null(fit_exp$loglik_mat))
    stop("method = \"elpd\" needs per-row log-likelihoods, which this model does not supply; ",
         "use method = \"taylor\" or method = \"cv\" with a `loss`", call. = FALSE)
  K <- length(fit_obs)
  p <- length(fit_exp$par)
  z <- matrix(stats::rnorm(n_draws * p), n_draws)
  score <- function(eta) {
    post <- mdpl_posterior(fit_exp, fit_obs, eta, robust = robust)
    Th <- sweep(z %*% chol(post$cov), 2, post$mean, "+")
    waic_elpd(fit_exp$loglik_mat(Th))
  }
  eta <- rep(0, K)
  tables <- vector("list", K)
  for (sweep_i in seq_len(if (K == 1L) 1L else 2L)) {
    for (k in seq_len(K)) {
      e <- vapply(grid, function(g) { x <- eta; x[k] <- g; score(x) }, numeric(1))
      eta[k] <- grid[which.max(e)]
      tables[[k]] <- data.frame(eta = grid, elpd = e)
    }
  }
  eta <- matrix(eta, K, 1L, dimnames = list(names(fit_obs) %||% paste0("source", seq_len(K)), "elpd"))
  structure(list(eta = eta, table = tables[[1]], tables = tables, K = K), class = "elpd_eta")
}

#' @export
print.elpd_eta <- function(x, digits = 3, ...) {
  cat("ELPD-selected borrowing weights (K =", x$K, "sources)\n")
  print(round(x$eta, digits))
  invisible(x)
}

# WAIC-based ELPD from an S x n matrix of log-likelihoods
waic_elpd <- function(ll) {
  m <- apply(ll, 2, max)
  lppd <- sum(m + log(colMeans(exp(sweep(ll, 2, m)))))
  lppd - sum(apply(ll, 2, stats::var))
}
