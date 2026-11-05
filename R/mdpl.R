#' Multi-source MDPL in one call
#'
#' Fits the outcome model on the trial and on each observational source,
#' selects the borrowing weights (one per source) and returns the
#' corresponding posteriors. This is the main entry point.
#'
#' @param model a model specification ([mdpl_glm()], [mdpl_causl()],
#'   [mdpl_custom()]) or a formula, which is shorthand for
#'   `mdpl_glm(formula, family)`.
#' @param data_exp trial data frame.
#' @param data_obs list of observational data frames (one weight each), or a
#'   single data frame (one weight). Names are carried to the output.
#' @param treat name of the numeric 0/1 treatment variable.
#' @param target optional extension, for `method = "taylor"` only. By default
#'   (`NULL`, or `"all"`) the borrowing is judged on the whole parameter vector,
#'   i.e. the full distribution. Set `"cate"` to judge it on the
#'   treatment-effect parameters only (nuisance parameters such as the
#'   intercept are then not penalised), or give explicit parameter names. The
#'   ELPD method always scores the whole trial likelihood, so `target` is
#'   ignored with a warning there.
#' @param method `"taylor"` (fast surrogate, [taylor_eta()]; returns one set of
#'   weights per trace multiplier), `"elpd"` (an in-sample ELPD proxy: WAIC on
#'   the trial data, [elpd_eta()]) or `"cv"` (gold-standard cross-validation on
#'   the trial, [cv_eta()]; options `type`, `folds`, `loss`, `grid` go through
#'   `...`). `"elpd"` needs per-row log-likelihoods, and `"cv"` needs them
#'   unless you give `loss` (see `per_row` in [mdpl_custom()]).
#' @param metric optional, `method = "taylor"` only: the matrix \eqn{H_0}
#'   (the Hessian of the validation risk) that scores the fused estimate, or a
#'   function of the trial model returning it. The default (`NULL`) is the
#'   trial loss, i.e. the trial information (or its `target` block). Overrides
#'   `target`.
#' @param trace_mult trace multipliers for `method = "taylor"`, default
#'   `c(1, 2)`. Both sets of weights are returned (columns `t1`, `t2` of
#'   `fit$eta`). `t2` (leave-one-out) is the default in every downstream
#'   function; `t1` (the more conservative Jensen bound) is kept as an option
#'   and is selected with `trace = "t1"`.
#' @param pooled if `TRUE`, stack the observational sources and use a single
#'   weight.
#' @param robust use the inverse sandwich covariance as the information.
#' @param bootstrap_CI if `TRUE`, also bootstrap the whole procedure (resample
#'   the trial and every source within treatment arm, refit, reselect the
#'   weights) `B` times and store the replicates. [mdpl_cate()], [mdpl_ate()]
#'   and [mdpl_coef()] then report bootstrap intervals by default. Default
#'   `FALSE`, because it multiplies the run time by about `B`.
#' @param B number of bootstrap replicates (at least 500 is advisable).
#' @param boot_seed optional seed for the bootstrap.
#' @param family used only when `model` is a formula.
#' @param cate `TRUE` (default) stores the CATE as a formula and function of the
#'   covariates in `estimates$cate` (see [mdpl_cate_formula()]; skipped quietly
#'   for models whose effect is not linear in the parameters). A data frame of
#'   covariate values (no treatment column needed) does that and also stores the
#'   CATE at those values in `estimates$cate_at`. `FALSE` or `NULL` stores
#'   nothing.
#' @param ate if `TRUE` (default), store the average treatment effect over the
#'   trial covariates in `estimates$ate`.
#' @param coef if `TRUE` (default), store the posterior summary of the
#'   treatment-effect coefficients in `estimates$coef`. All MDPL models are
#'   parametric, so these exist for every model; set `FALSE` to leave them out
#'   (for example when the parameters have no direct meaning).
#' @param ... passed to [elpd_eta()] (e.g. `grid`, `n_draws`) or [cv_eta()].
#'
#' @return an object of class `"mdpl"` with `eta` (matrix, one column per
#'   trace), `posterior` (list of [mdpl_posterior()] named like the columns of
#'   `eta`), `rct` (the trial-only posterior), `fit_exp`, `fit_obs`,
#'   `selection`, `estimates` (a list with `trace` and whichever of `coef`,
#'   `ate`, `cate`, `cate_at` were requested; `NULL` if none) and the inputs. With `bootstrap_CI = TRUE`,
#'   also `boot`, the bootstrap replicates of the fitted parameters.
#' @examples
#' set.seed(2)
#' sim <- function(n, bias = 0) {
#'   X <- rnorm(n, 1); U <- rnorm(n)
#'   A <- rbinom(n, 1, if (bias == 0) 0.5 else plogis(X - 1 + bias * U))
#'   data.frame(X = X, A = A, Y = 0.5 + X + 0.2 * A + 0.1 * A * X + bias * U + rnorm(n))
#' }
#' fit <- mdpl(Y ~ A * X, sim(300), lapply(c(0, 1), function(b) sim(1000, b)))
#' fit
#' fit$eta                      # selected weights
#' fit$estimates$coef           # treatment-effect coefficients
#' fit$estimates$ate            # average treatment effect
#' fit2 <- mdpl(Y ~ A * X, sim(300), lapply(c(0, 1), function(b) sim(1000, b)),
#'              cate = data.frame(X = c(0, 1, 2)))
#' fit2$estimates$cate_at
#' fit$estimates$cate                                          # formula, callable
#' fit$estimates$cate(data.frame(X = c(0, 1, 2)))
#' mdpl_cate(fit, data.frame(X = c(0, 1, 2)))                 # uses t2
#' mdpl_cate(fit, data.frame(X = c(0, 1, 2)), trace = "t1")   # more conservative
#' @export
mdpl <- function(model, data_exp, data_obs, treat = "A", target = NULL,
                 method = c("taylor", "elpd", "cv"), trace_mult = c(1, 2),
                 pooled = FALSE, robust = FALSE, bootstrap_CI = FALSE, B = 500L,
                 boot_seed = NULL, family = stats::gaussian(), metric = NULL,
                 cate = TRUE, ate = TRUE, coef = TRUE, ...) {
  method <- match.arg(method)
  if (method != "taylor" && !is.null(target) && !identical(target, "all"))
    warning("`target` only applies to method = \"taylor\"; the ", method, " method scores the trial loss")
  if (method != "taylor" && !is.null(metric))
    warning("`metric` only applies to method = \"taylor\"")
  spec <- as_spec(model, family)
  if (is.data.frame(data_obs)) data_obs <- list(data_obs)
  res <- run_mdpl(spec, as.data.frame(data_exp), lapply(data_obs, as.data.frame), treat, target,
                  method, trace_mult, pooled, robust, metric, ...)
  res$data_exp <- data_exp
  res$data_obs <- data_obs
  res$dots <- list(...)
  if (bootstrap_CI) {
    res$boot <- boot_theta(res, B, NULL, boot_seed)
    res$boot$B <- B
  }
  res$estimates <- collect_estimates(res, cate, ate, coef)
  res
}

run_mdpl <- function(spec, data_exp, data_obs, treat, target, method, trace_mult,
                     pooled, robust, metric = NULL, ...) {
  fit_exp <- spec$fit(data_exp)
  fit_obs <- if (pooled) list(spec$fit(do.call(rbind, data_obs)))
             else lapply(data_obs, spec$fit)
  if (!is.null(names(data_obs)) && !pooled) names(fit_obs) <- names(data_obs)
  check_same_target(fit_exp, fit_obs)
  tg <- if (identical(target, "cate")) treatment_terms(fit_exp, treat)
        else if (identical(target, "all")) NULL else target
  A <- if (is.null(metric)) NULL else {
    H <- if (is.function(metric)) metric(fit_exp) else metric
    if (!is.matrix(H) || any(dim(H) != length(fit_exp$par)))
      stop("`metric` must be a p x p matrix (or a function returning one)")
    -(H + t(H)) / 2
  }
  sel <- switch(method,
    taylor = taylor_eta(fit_exp, fit_obs, target = tg, A = A, trace_mult = trace_mult, robust = robust),
    elpd = elpd_eta(fit_exp, fit_obs, robust = robust, ...),
    cv = cv_eta(fit_exp, fit_obs, spec = spec, data_exp = data_exp, robust = robust, ...))
  post <- stats::setNames(lapply(seq_len(ncol(sel$eta)), function(j)
    mdpl_posterior(fit_exp, fit_obs, sel$eta[, j], treat = treat, robust = robust)),
    colnames(sel$eta))
  structure(list(eta = sel$eta, posterior = post,
                 rct = mdpl_posterior(fit_exp, fit_obs, rep(0, length(fit_obs)),
                                      treat = treat, robust = robust),
                 fit_exp = fit_exp, fit_obs = fit_obs, selection = sel,
                 spec = spec, treat = treat, target = if (method == "taylor") tg else NULL, method = method,
                 trace_mult = trace_mult, pooled = pooled, robust = robust,
                 metric = metric),
            class = "mdpl")
}

collect_estimates <- function(fit, cate, ate, coef) {
  want_cate <- is.data.frame(cate) || isTRUE(cate)
  if (!want_cate && !ate && !coef) return(NULL)
  tr <- default_trace(fit)
  get <- function(what, off, expr) tryCatch(expr, error = function(e) {
    warning("Could not compute the ", what, ": ", conditionMessage(e),
            "\nTurn it off with `", off, "`.", call. = FALSE)
    NULL
  })
  out <- list(trace = tr)
  if (coef) out$coef <- get("treatment-effect coefficients", "coef = FALSE",
                            mdpl_coef(fit, trace = tr))
  if (ate) out$ate <- get("ATE", "ate = FALSE", mdpl_ate(fit, trace = tr))
  if (want_cate) {
    out$cate <- tryCatch(mdpl_cate_formula(fit, trace = tr), error = function(e) NULL)
    if (is.data.frame(cate))
      out$cate_at <- get("CATE", "cate = TRUE", mdpl_cate(fit, cate, trace = tr))
  }
  out
}

#' @export
print.mdpl <- function(x, digits = 3, ...) {
  cat("MDPL fit (", x$spec$label, ")\n", sep = "")
  cat(length(x$fit_obs), " observational source(s); weights chosen by ", x$method,
      "; scored on ", score_desc(x), "\n", sep = "")
  cat("Borrowing weights eta",
      if (all(c("t1", "t2") %in% colnames(x$eta))) " (t2 is the default; t1 is more conservative)" else "",
      ":\n", sep = "")
  print(round(x$eta, digits))
  if (!is.null(x$estimates)) {
    e <- x$estimates
    if (!is.null(e$coef)) {
      cat("Treatment-effect coefficients (", e$trace, "):\n", sep = "")
      print(format(e$coef, digits = digits), row.names = FALSE)
    }
    if (!is.null(e$ate)) { cat("ATE:\n"); print(format(e$ate, digits = digits), row.names = FALSE) }
    if (!is.null(e$cate)) cat(attr(e$cate, "text"), "\n", sep = "")
    if (!is.null(e$cate_at)) {
      cat("CATE at the requested values:\n")
      print(format(e$cate_at, digits = digits), row.names = FALSE)
    }
  }
  if (!is.null(x$boot)) cat("Bootstrap intervals stored (B = ", x$boot$B, ")\n", sep = "")
  invisible(x)
}

#' Description of what the weights are scored on
#' @param x an `mdpl` fit.
#' @return a string.
#' @keywords internal
#' @export
score_desc <- function(x) {
  if (x$method == "elpd") return("the trial log-likelihood (ELPD proxy)")
  if (x$method == "cv") return("held-out trial units")
  if (!is.null(x$metric)) return("a user-supplied metric")
  if (is.null(x$target)) "all parameters" else paste(x$target, collapse = ", ")
}
