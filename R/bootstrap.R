# B x p x traces array of posterior means from bootstrapping the whole procedure
boot_theta <- function(x, B, trace, seed) {
  if (!inherits(x, "mdpl")) stop("`x` must be an `mdpl` object")
  if (!is.null(seed)) set.seed(seed)
  traces <- trace %||% names(x$posterior)
  if (!all(traces %in% names(x$posterior))) stop("unknown `trace`")
  tm <- x$trace_mult[match(traces, names(x$posterior))]
  treat <- x$treat
  resample <- function(d) {
    d <- as.data.frame(d)
    idx <- unlist(lapply(split(seq_len(nrow(d)), d[[treat]]),
                         function(g) g[sample.int(length(g), replace = TRUE)]),
                  use.names = FALSE)
    d[idx, , drop = FALSE]
  }
  p <- length(x$posterior[[1]]$mean)
  draws <- array(NA_real_, c(B, p, length(traces)),
                 dimnames = list(NULL, names(x$posterior[[1]]$mean), traces))
  for (b in seq_len(B)) {
    fb <- tryCatch(
      do.call(run_mdpl, c(list(x$spec, resample(x$data_exp), lapply(x$data_obs, resample),
                               treat, x$target, x$method, tm, x$pooled, x$robust, x$metric), x$dots)),
      error = function(e) NULL)
    if (is.null(fb) || ncol(fb$eta) != length(traces)) next
    for (j in seq_along(traces)) draws[b, , j] <- fb$posterior[[j]]$mean
  }
  list(draws = draws, traces = traces)
}

# stored replicates if B is not given, otherwise a fresh bootstrap
get_boot <- function(x, B, trace, seed) {
  if (is.null(B) && !is.null(x$boot)) {
    traces <- trace %||% default_trace(x)
    if (!all(traces %in% names(x$posterior))) stop("unknown `trace`")
    return(list(draws = x$boot$draws[, , traces, drop = FALSE], traces = traces))
  }
  boot_theta(x, B %||% 500L, trace %||% default_trace(x), seed)
}

# interval summaries; `fun(Theta)` maps S x p parameter vectors to S x m effects
finish_boot <- function(x, bt, fun, level, interval, labels) {
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  B <- dim(bt$draws)[1]
  m <- ncol(fun(x$posterior[[1]]$mean))
  dr <- array(NA_real_, c(B, m, length(bt$traces)), dimnames = list(NULL, NULL, bt$traces))
  out <- do.call(rbind, lapply(seq_along(bt$traces), function(j) {
    ok <- stats::complete.cases(bt$draws[, , j])
    d <- matrix(NA_real_, B, m)
    d[ok, ] <- fun(matrix(bt$draws[ok, , j], ncol = dim(bt$draws)[2]))
    dr[, , j] <<- d
    est <- fun(x$posterior[[bt$traces[j]]]$mean)[1, ]
    se <- apply(d, 2, stats::sd, na.rm = TRUE)
    if (interval == "percentile") {
      lo <- apply(d, 2, stats::quantile, a, na.rm = TRUE, names = FALSE)
      hi <- apply(d, 2, stats::quantile, 1 - a, na.rm = TRUE, names = FALSE)
    } else { lo <- est - z * se; hi <- est + z * se }
    res <- data.frame(trace = bt$traces[j], estimate = est, se = se, lower = lo, upper = hi)
    if (is.null(labels)) res else cbind(labels, res)
  }))
  rownames(out) <- NULL
  attr(out, "draws") <- dr
  out
}

#' Bootstrap intervals for the CATE, reselecting the weights each time
#'
#' The posterior interval from [mdpl_cate()] conditions on the selected
#' weights, so it ignores the uncertainty in `eta` and can under-cover when
#' the observational sources are biased. This bootstraps the whole procedure:
#' in each replicate the trial and every source are resampled within treatment
#' arm, the models are refitted, the weights are reselected and the effect is
#' recomputed. With `method = "taylor"` a few hundred replicates take seconds
#' for a GLM.
#'
#' @param x an `mdpl` object.
#' @param newdata covariate values at which to evaluate the CATE.
#' @param B number of bootstrap replicates (use at least 500 for intervals). If
#'   `NULL` (default), the replicates stored by `mdpl(..., bootstrap_CI = TRUE)`
#'   are reused when present, otherwise 500 new ones are drawn.
#' @param scale,level as in [mdpl_cate()].
#' @param trace which trace multiplier: `"t2"` (default) or `"t1"`.
#' @param seed optional seed.
#' @param interval `"percentile"` or `"normal"` (`estimate +/- z * se`).
#' @return a data frame with one row per trace and row of `newdata`: the
#'   covariates, `trace`, `estimate` (from the original data), `se` (bootstrap
#'   standard deviation), `lower`, `upper`. The replicates are in
#'   `attr(, "draws")` (`B x nrow(newdata) x traces`).
#' @examples
#' \donttest{
#' set.seed(3)
#' sim <- function(n, bias = 0) {
#'   X <- rnorm(n, 1); U <- rnorm(n)
#'   A <- rbinom(n, 1, if (bias == 0) 0.5 else plogis(X - 1 + bias * U))
#'   data.frame(X = X, A = A, Y = 0.5 + X + 0.2 * A + 0.1 * A * X + bias * U + rnorm(n))
#' }
#' fit <- mdpl(Y ~ A * X, sim(300), lapply(c(0, 1), function(b) sim(1000, b)))
#' mdpl_cate_boot(fit, data.frame(X = c(0, 1, 2)), B = 50, trace = "t2")
#' }
#' @export
mdpl_cate_boot <- function(x, newdata, B = NULL, scale = c("response", "link"),
                           level = 0.95, trace = NULL, seed = NULL,
                           interval = c("percentile", "normal")) {
  scale <- match.arg(scale); interval <- match.arg(interval)
  bt <- get_boot(x, B, trace, seed)
  nd <- prep_newdata(newdata, x$treat)
  fun <- function(Th) effect_matrix(x$fit_exp, Th, nd, x$treat, scale, FALSE)
  labels <- nd[, setdiff(names(nd), x$treat), drop = FALSE]
  finish_boot(x, bt, fun, level, interval, labels)
}

#' Bootstrap interval for the average treatment effect
#'
#' Same procedure as [mdpl_cate_boot()] for the ATE over the covariate
#' distribution of `data` (default: the trial), which is held fixed in every
#' replicate.
#'
#' @inheritParams mdpl_cate_boot
#' @param data data frame defining the target population.
#' @return a data frame with one row per trace: `trace`, `estimate`, `se`,
#'   `lower`, `upper`; the replicates are in `attr(, "draws")`.
#' @export
mdpl_ate_boot <- function(x, data = x$data_exp, B = NULL, scale = c("response", "link"),
                          level = 0.95, trace = NULL, seed = NULL,
                          interval = c("percentile", "normal")) {
  scale <- match.arg(scale); interval <- match.arg(interval)
  bt <- get_boot(x, B, trace, seed)
  nd <- prep_newdata(data, x$treat)
  fun <- function(Th) effect_matrix(x$fit_exp, Th, nd, x$treat, scale, TRUE)
  finish_boot(x, bt, fun, level, interval, NULL)
}

#' Bootstrap intervals for the treatment-effect parameters
#'
#' Bootstrap (as in [mdpl_cate_boot()]) for the parameters that make up the
#' effect, e.g. `A` and `A:X`.
#'
#' @inheritParams mdpl_cate_boot
#' @param terms parameter names; default the treatment parameters.
#' @param level interval level.
#' @return a data frame with one row per trace and parameter: `term`, `trace`,
#'   `estimate`, `se`, `lower`, `upper`. `attr(, "draws")` holds the replicates.
#' @export
mdpl_coef_boot <- function(x, terms = NULL, B = NULL, level = 0.95, trace = NULL,
                           seed = NULL, interval = c("percentile", "normal")) {
  interval <- match.arg(interval)
  bt <- get_boot(x, B, trace, seed)
  terms <- terms %||% treatment_terms(x$fit_exp, x$treat)
  ix <- match(terms, names(x$fit_exp$par))
  fun <- function(Th) matrix(Th, ncol = length(x$fit_exp$par))[, ix, drop = FALSE]
  out <- finish_boot(x, bt, fun, level, interval, NULL)
  dr <- attr(out, "draws")
  out <- cbind(term = rep(terms, times = nrow(out) / length(terms)), out)
  attr(out, "draws") <- dr
  out
}
