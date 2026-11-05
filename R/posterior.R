#' MDPL power-likelihood posterior
#'
#' The power posterior of the parameters, in its normal approximation, is
#' \eqn{N(m_\eta, S_\eta^{-1})} with \eqn{S_\eta = I_e + \sum_k \eta_k I_k}
#' and \eqn{m_\eta = S_\eta^{-1}(I_e\hat\theta_e + \sum_k \eta_k I_k
#' \hat\theta_k)}, where \eqn{I_e, I_k} are the information matrices of the
#' trial and of source \eqn{k}. For the Gaussian linear model this is exact.
#'
#' @param fit_exp trial [mdpl_model][new_mdpl_model].
#' @param fit_obs one model or a list of them (one per source).
#' @param eta weights in `[0, 1]`, one per source.
#' @param treat name of the treatment variable (stored for the effect
#'   functions).
#' @param robust use the inverse sandwich covariance instead of the
#'   model-based information.
#' @return an object of class `"mdpl_posterior"` with `mean`, `cov`, `eta`,
#'   `model` (the trial model) and `treat`.
#' @export
mdpl_posterior <- function(fit_exp, fit_obs, eta, treat = NULL, robust = FALSE) {
  if (inherits(fit_obs, "mdpl_model")) fit_obs <- list(fit_obs)
  check_same_target(fit_exp, fit_obs)
  eta <- as.numeric(eta)
  if (length(eta) != length(fit_obs)) stop("`eta` needs one weight per observational fit")
  eta <- pmin(pmax(eta, 0), 1)
  Ie <- model_info(fit_exp, robust)
  S <- Ie
  b <- Ie %*% fit_exp$par
  for (k in seq_along(fit_obs)) {
    Ik <- model_info(fit_obs[[k]], robust)
    S <- S + eta[k] * Ik
    b <- b + eta[k] * (Ik %*% fit_obs[[k]]$par)
  }
  V <- solve(S)
  V <- (V + t(V)) / 2
  nm <- names(fit_exp$par)
  structure(list(mean = stats::setNames(as.numeric(V %*% b), nm),
                 cov = `dimnames<-`(V, list(nm, nm)),
                 eta = eta, model = fit_exp, treat = treat),
            class = "mdpl_posterior")
}

#' @export
print.mdpl_posterior <- function(x, ...) {
  cat("MDPL posterior; eta =", paste(round(x$eta, 3), collapse = ", "), "\n")
  print(signif(x$mean, 4))
  invisible(x)
}

#' Draw from an MDPL posterior
#'
#' @param post an `mdpl_posterior`.
#' @param n number of draws.
#' @return an `n x p` matrix, one draw per row.
#' @export
mdpl_draws <- function(post, n = 1000L) {
  R <- chol(post$cov)
  z <- matrix(stats::rnorm(n * length(post$mean)), n)
  sweep(z %*% R, 2, post$mean, "+")
}
