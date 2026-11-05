#' Specification: generalised linear model
#'
#' Supported families are `gaussian()`, `binomial()` (logit, probit, cloglog),
#' `poisson()` and `Gamma()`. The same formula and family are used for the
#' trial and every observational source. An `offset()` term is supported (for
#' example `offset(log(followup))` in a Poisson rate model); a response-scale
#' effect then needs the offset variable in `newdata`. The treatment must be a numeric 0/1
#' variable on the right-hand side; interactions with covariates give a
#' conditional effect.
#'
#' The information matrix is the expected information
#' \eqn{X^\top W X / \phi} at the fitted values, with dispersion \eqn{\phi}
#' (the ML residual variance for the Gaussian family, the Pearson estimate for
#' Gamma, and 1 for binomial and Poisson). The robust covariance is the usual
#' sandwich. For the Gaussian family the likelihood used for the ELPD fixes
#' \eqn{\phi} at its estimate.
#'
#' @param formula outcome model, e.g. `Y ~ A * X`.
#' @param family a family object or its name.
#' @return an object of class `"mdpl_spec"`.
#' @examples
#' spec <- mdpl_glm(y ~ A * x, binomial())
#' spec
#' @export
mdpl_glm <- function(formula, family = stats::gaussian()) {
  if (is.character(family)) family <- get(family, mode = "function")()
  if (is.function(family)) family <- family()
  if (!family$family %in% c("gaussian", "binomial", "poisson", "Gamma"))
    stop("family must be gaussian, binomial, poisson or Gamma; use mdpl_custom() for others")
  structure(list(fit = function(data) fit_glm_model(formula, family, data),
                 formula = formula, family = family,
                 label = paste0("GLM, ", family$family, "(", family$link, "): ",
                                paste(deparse(formula), collapse = ""))),
            class = "mdpl_spec")
}

fit_glm_model <- function(formula, family, data) {
  tt <- stats::terms(formula)
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  if (is.factor(y)) y <- as.numeric(y != levels(y)[1])
  if (is.logical(y)) y <- as.numeric(y)
  X <- stats::model.matrix(tt, mf)
  off <- stats::model.offset(mf) %||% rep(0, nrow(X))
  n <- nrow(X); p <- ncol(X)
  if (n <= p) stop("Need more rows than coefficients")
  if (family$family == "binomial" && !all(y %in% c(0, 1))) stop("binomial outcome must be 0/1")
  fit <- suppressWarnings(stats::glm.fit(X, y, family = family, offset = off,
                                         control = stats::glm.control(epsilon = 1e-10, maxit = 100)))
  if (fit$rank < p) stop("Outcome design matrix is rank deficient")
  if (!fit$converged) warning("glm.fit did not converge")
  mu <- fit$fitted.values; eta <- fit$linear.predictors
  v <- family$variance(mu)
  phi <- switch(family$family,
                gaussian = sum((y - mu)^2) / n,
                Gamma = sum((y - mu)^2 / v) / (n - p),
                1)
  if (!is.finite(phi) || phi <= 0) stop("Dispersion is not positive")
  Xw <- X * sqrt(fit$weights)
  OI <- crossprod(Xw) / phi
  score <- X * ((y - mu) * family$mu.eta(eta) / (phi * v))
  H <- solve(OI)
  sand <- H %*% crossprod(score) %*% H
  par <- stats::setNames(fit$coefficients, colnames(X))
  tt_rhs <- stats::delete.response(tt)
  xlev <- stats::.getXlevels(tt, mf)

  design <- function(newdata, treat, value) {
    nd <- as.data.frame(newdata)
    nd[[treat]] <- value
    mfn <- stats::model.frame(tt_rhs, nd, xlev = xlev)
    M <- stats::model.matrix(tt_rhs, mfn)
    attr(M, "offset") <- stats::model.offset(mfn) %||% rep(0, nrow(M))
    if (!identical(colnames(M), names(par)))
      stop("`newdata` does not reproduce the coefficients of the fitted model")
    M
  }
  linkinv <- family$linkinv
  ll_one <- switch(family$family,
    gaussian = function(mu) stats::dnorm(y, mu, sqrt(phi), log = TRUE),
    binomial = function(mu) {
      mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
      y * log(mu) + (1 - y) * log1p(-mu)
    },
    poisson = function(mu) stats::dpois(y, mu, log = TRUE),
    Gamma = function(mu) stats::dgamma(y, shape = 1 / phi, scale = mu * phi, log = TRUE))
  loglik_mat <- function(Theta) {
    Theta <- matrix(Theta, ncol = p)
    t(vapply(seq_len(nrow(Theta)), function(s) ll_one(linkinv(X %*% Theta[s, ] + off)),
             numeric(n)))
  }
  effect <- function(Theta, newdata, treat, scale = "response") {
    Theta <- matrix(Theta, ncol = p)
    d1 <- design(newdata, treat, 1); d0 <- design(newdata, treat, 0)
    e1 <- sweep(d1 %*% t(Theta), 1, attr(d1, "offset"), "+")
    e0 <- sweep(d0 %*% t(Theta), 1, attr(d0, "offset"), "+")
    t(if (scale == "link") e1 - e0 else linkinv(e1) - linkinv(e0))
  }
  contrast <- function(newdata, treat, scale = "response") {
    if (scale == "link" || family$link == "identity")
      { L <- design(newdata, treat, 1) - design(newdata, treat, 0); attr(L, "offset") <- NULL; L }
    else NULL
  }
  new_mdpl_model(par, OI, sand, n, loglik_mat, effect, contrast,
                 treatment_terms = function(treat) {
                   hit <- vapply(strsplit(names(par), ":", fixed = TRUE),
                                 function(z) treat %in% z, logical(1))
                   names(par)[hit]
                 },
                 family = family, formula = formula, phi = phi,
                 linear = if (family$family == "gaussian" && family$link == "identity")
                   list(G = X, y = y - off, sigma2 = phi),
                 label = paste0("glm/", family$family))
}
