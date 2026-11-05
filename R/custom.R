#' Specification: user-supplied likelihood
#'
#' For any smooth parametric model not covered by [mdpl_glm()] or
#' [mdpl_causl()]. You supply the per-row log-likelihood and a starting
#' value; the package finds the maximum likelihood estimate with `optim`, and
#' gets the information matrix and the sandwich covariance from numerical
#' derivatives (package numDeriv), so no gradients are needed. Run
#' [mdpl_check_model()] first: it tests the conditions the fusion relies on.
#'
#' The model must be a regular parametric likelihood: finite-dimensional
#' parameter in the interior of its space, identifiable, a log-likelihood that
#' is smooth in the parameter, and the same parameterisation for the trial and
#' every observational source.
#'
#' @param loglik `function(theta, data)` returning the vector of per-row
#'   log-likelihoods (the row is the unit; for clustered data make each row a
#'   cluster). With `per_row = FALSE` it returns the total log-likelihood
#'   (a single number) instead.
#' @param per_row `TRUE` (default) if `loglik` returns per-row values. If
#'   `FALSE`, `loglik` returns the total log-likelihood. The model is then
#'   fitted and fused as usual (information from the numerical Hessian), but
#'   there are no per-row values, so there is no sandwich covariance
#'   (`robust = TRUE` is unavailable), and `method = "elpd"` and the default
#'   validation loss of `method = "cv"` cannot be used; give `loss` to
#'   `method = "cv"`. `method = "taylor"` works as normal.
#' @param start numeric starting value (named, to get named parameters).
#' @param effect `function(theta, newdata, treat)` returning the treatment
#'   effect at each row of `newdata` for ONE parameter vector `theta`. Needed
#'   for CATE and ATE. If the effect is linear in `theta` supply `contrast`
#'   instead (or as well) for exact intervals.
#' @param contrast `function(newdata, treat)` returning a matrix `L` with
#'   effect `= L %*% theta`.
#' @param treatment_terms `function(treat)` giving the names of the
#'   parameters that determine the effect; default all parameters.
#' @param optim_method method passed to [stats::optim()].
#' @param label text used when printing.
#' @return an object of class `"mdpl_spec"`.
#' @examples
#' # Poisson regression written by hand
#' ll <- function(theta, data) stats::dpois(data$y, exp(theta[1] + theta[2] * data$A), log = TRUE)
#' spec <- mdpl_custom(ll, start = c(a = 0, A = 0),
#'                     contrast = function(newdata, treat) matrix(c(0, 1), nrow(newdata), 2, byrow = TRUE))
#' set.seed(1)
#' d <- data.frame(A = rbinom(200, 1, .5)); d$y <- rpois(200, exp(0.2 + 0.5 * d$A))
#' mdpl_fit(spec, d)
#' @export
mdpl_custom <- function(loglik, start, effect = NULL, contrast = NULL,
                        treatment_terms = NULL, optim_method = "BFGS",
                        label = "custom likelihood", per_row = TRUE) {
  structure(list(
    fit = function(data) fit_custom_model(loglik, start, effect, contrast,
                                          treatment_terms, optim_method, data, per_row),
    label = label, loglik = loglik, start = start), class = "mdpl_spec")
}

fit_custom_model <- function(loglik, start, effect, contrast, treatment_terms, method, data, per_row = TRUE) {
  p <- length(start)
  nm <- names(start) %||% paste0("theta", seq_len(p))
  f <- function(th) { v <- loglik(th, data); if (any(!is.finite(v))) -1e10 else sum(v) }
  opt <- stats::optim(start, f, method = method, control = list(fnscale = -1, maxit = 1000, reltol = 1e-12))
  if (opt$convergence != 0) warning("optim did not converge")
  th <- stats::setNames(opt$par, nm)
  rows <- function(t) as.numeric(loglik(t, data))
  if (per_row) {
    n <- length(rows(th))
    Hs <- -numDeriv::hessian(function(t) sum(rows(t)), th)
    J <- numDeriv::jacobian(rows, th)
  } else {
    if (length(rows(th)) != 1L)
      stop("with `per_row = FALSE`, `loglik` must return a single total log-likelihood")
    n <- nrow(data)
    Hs <- -numDeriv::hessian(rows, th)
    J <- NULL
  }
  Hs <- (Hs + t(Hs)) / 2
  sand <- if (is.null(J)) NULL else
    tryCatch({ H <- solve(Hs); H %*% crossprod(J) %*% H }, error = function(e) NULL)

  loglik_mat <- if (!per_row) NULL else function(Theta) {
    Theta <- matrix(Theta, ncol = p)
    t(vapply(seq_len(nrow(Theta)), function(s) rows(Theta[s, ]), numeric(n)))
  }
  eff <- if (is.null(effect)) function(Theta, newdata, treat, scale = "response")
    stop("this model has no `effect` function; supply `effect` or `contrast` to mdpl_custom()")
  else function(Theta, newdata, treat, scale = "response") {
    Theta <- matrix(Theta, ncol = p)
    t(vapply(seq_len(nrow(Theta)), function(s) as.numeric(effect(Theta[s, ], newdata, treat)),
             numeric(nrow(newdata))))
  }
  con <- if (is.null(contrast)) function(newdata, treat, scale = "response") NULL
         else function(newdata, treat, scale = "response") contrast(newdata, treat)
  if (is.null(effect) && !is.null(contrast))
    eff <- function(Theta, newdata, treat, scale = "response")
      matrix(Theta, ncol = p) %*% t(contrast(newdata, treat))
  new_mdpl_model(th, Hs, sand, n, loglik_mat, eff, con,
                 treatment_terms = treatment_terms %||% function(treat) names(th),
                 label = "custom")
}

#' Check that a likelihood is suitable for MDPL
#'
#' Fits the specification to `data` and tests the conditions the fusion
#' machinery relies on: finite per-row log-likelihood, a vanishing gradient at
#' the estimate, a positive-definite information matrix with a reasonable
#' condition number, a score that matches the numerical Hessian (information
#' identity, a warning sign for misspecification), and parameters away from the
#' boundary of their support. Pass the trial and an observational data set
#' together to also check that the parameter vectors agree in length.
#'
#' @param spec a specification (usually [mdpl_custom()]).
#' @param data a data frame, or a list of data frames (e.g. trial and source).
#' @param stop_on_fail signal an error if any check fails.
#' @return a data frame of checks with columns `check`, `passed`, `detail`
#'   (invisibly, after printing).
#' @export
mdpl_check_model <- function(spec, data, stop_on_fail = FALSE) {
  spec <- as_spec(spec)
  data <- if (is.data.frame(data)) list(data) else data
  rows <- list()
  add <- function(check, passed, detail = "") rows[[length(rows) + 1L]] <<-
    data.frame(check = check, passed = isTRUE(passed), detail = detail)
  fits <- lapply(data, function(d) tryCatch(spec$fit(as.data.frame(d)), error = function(e) e))
  bad <- vapply(fits, inherits, logical(1), "error")
  add("model fits", !any(bad), if (any(bad)) conditionMessage(fits[[which(bad)[1]]]) else "")
  if (!any(bad)) {
    m <- fits[[1]]
    if (is.null(m$loglik_mat)) {
      add("per-row checks skipped", TRUE,
          "no per-row log-likelihood: ELPD weights, the sandwich and the default CV loss are unavailable")
      ll <- NULL
    } else ll <- m$loglik_mat(matrix(m$par, 1))
    if (!is.null(ll)) add("per-row log-likelihood is finite", all(is.finite(ll)), paste(ncol(ll), "rows"))
    if (!is.null(ll) && all(is.finite(ll))) {
      g <- tryCatch(numDeriv::jacobian(function(t) as.numeric(m$loglik_mat(matrix(t, 1))), m$par),
                    error = function(e) NULL)
      if (!is.null(g)) {
        gs <- colSums(g)
        add("gradient is zero at the estimate", max(abs(gs)) < 1e-2 * sqrt(m$n),
            paste("max |score| =", signif(max(abs(gs)), 3)))
      }
    }
    ev <- eigen(m$OI, symmetric = TRUE, only.values = TRUE)$values
    add("information is positive definite", all(ev > 0), paste("min eigenvalue", signif(min(ev), 3)))
    add("information is well conditioned", all(ev > 0) && max(ev) / min(ev) < 1e10,
        paste("condition number", signif(max(ev) / max(min(ev), 1e-300), 3)))
    if (!is.null(m$sandwich) && all(ev > 0)) {
      r <- sqrt(diag(m$sandwich)) / sqrt(diag(solve(m$OI)))
      add("model-based and robust standard errors agree", all(r > 0.5 & r < 2),
          paste("ratio range", paste(signif(range(r), 3), collapse = " - ")))
    }
    add("parameters have the same length in every data set",
        length(unique(vapply(fits, function(f) length(f$par), 1L))) == 1L)
    if (length(fits) > 1L)
      add("parameter names agree", length(unique(lapply(fits, function(f) names(f$par)))) == 1L)
  }
  out <- do.call(rbind, rows)
  print(out, row.names = FALSE)
  if (stop_on_fail && !all(out$passed)) stop("model failed the checks above")
  invisible(out)
}
