#' MDPL model objects
#'
#' Every model family supported by MDPL is reduced to the same small object,
#' an `mdpl_model`, which is all the fusion machinery sees. Users normally do
#' not build these by hand: create a model *specification* with [mdpl_glm()],
#' [mdpl_causl()] or [mdpl_custom()] and fit it with [mdpl_fit()].
#'
#' An `mdpl_model` is a list with
#' \describe{
#'   \item{`par`}{named parameter estimate \eqn{\hat\theta}.}
#'   \item{`OI`}{information matrix (observed or expected) at `par`.}
#'   \item{`sandwich`}{robust covariance of `par`, or `NULL`.}
#'   \item{`robust_info`}{optional; if present it is used as the robust
#'     information instead of `solve(sandwich)` (for models whose information
#'     is singular because some parameters carry none).}
#'   \item{`n`}{number of units.}
#'   \item{`loglik_mat(Theta)`}{`Theta` is an `S x p` matrix of parameter
#'     vectors; returns the `S x n` matrix of per-unit log-likelihoods. May be
#'     `NULL` for a model without per-row log-likelihoods; such a model can
#'     be fused with `method = "taylor"`, or `"cv"` with a custom `loss`, but
#'     not with `method = "elpd"`.}
#'   \item{`effect(Theta, newdata, treat, scale)`}{treatment effect at each row
#'     of `newdata` for each row of `Theta` (an `S x nrow(newdata)` matrix).}
#'   \item{`contrast(newdata, treat, scale)`}{matrix `L` with effect `= L %*%
#'     theta` if the effect is linear in `theta`, otherwise `NULL`.}
#'   \item{`treatment_terms(treat)`}{names of the parameters that determine the
#'     treatment effect.}
#' }
#'
#' @param par,OI,sandwich,n,loglik_mat,effect,contrast,treatment_terms see
#'   Details.
#' @param ... further elements stored in the object (e.g. the family).
#' @return an object of class `"mdpl_model"`.
#' @name mdpl_model
#' @export
new_mdpl_model <- function(par, OI, sandwich = NULL, n, loglik_mat, effect,
                           contrast = function(newdata, treat, scale) NULL,
                           treatment_terms = function(treat) names(par), ...) {
  p <- length(par)
  if (is.null(names(par))) names(par) <- paste0("theta", seq_len(p))
  if (!is.matrix(OI) || any(dim(OI) != p)) stop("`OI` must be a p x p matrix")
  dimnames(OI) <- list(names(par), names(par))
  if (!is.null(sandwich)) dimnames(sandwich) <- list(names(par), names(par))
  structure(list(par = par, OI = OI, sandwich = sandwich, n = n,
                 loglik_mat = loglik_mat, effect = effect, contrast = contrast,
                 treatment_terms = treatment_terms, ...),
            class = "mdpl_model")
}

#' @export
print.mdpl_model <- function(x, ...) {
  cat("<mdpl_model> ", x$label %||% "", " n = ", x$n, ", ", length(x$par), " parameters\n", sep = "")
  print(signif(x$par, 4))
  invisible(x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# comparable summary of a model family: family objects carry closures, so use their names
family_key <- function(x) {
  if (inherits(x, "family")) c(x$family, x$link)
  else if (is.list(x)) lapply(x, family_key)
  else if (is.function(x)) NULL
  else x
}

#' Check that a source estimates the same thing as the trial
#'
#' Sources can only be merged with the trial if they target the same
#' parameter vector: the same parameter names in the same order and, where the
#' models record one, the same model family (and estimand, such as the event
#' type of a competing-risks model). Anything else is an error.
#'
#' @param fit_exp trial model.
#' @param fit_obs one model or a list of them.
#' @return `TRUE` invisibly.
#' @keywords internal
#' @export
check_same_target <- function(fit_exp, fit_obs) {
  if (inherits(fit_obs, "mdpl_model")) fit_obs <- list(fit_obs)
  for (k in seq_along(fit_obs)) {
    f <- fit_obs[[k]]
    if (!identical(names(f$par), names(fit_exp$par)))
      stop("Source ", k, " does not target the same parameters as the trial (trial: ",
           paste(names(fit_exp$par), collapse = ", "), "; source: ",
           paste(names(f$par), collapse = ", "), "). Studies that target different ",
           "quantities cannot be merged.", call. = FALSE)
    if (!identical(family_key(f$family), family_key(fit_exp$family)))
      stop("Source ", k, " uses a different model family or estimand from the trial; ",
           "studies that target different quantities cannot be merged.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Information used for fusion
#'
#' Model-based information, or the inverse robust covariance.
#' @param m a model from [new_mdpl_model()].
#' @param robust use the inverse sandwich covariance.
#' @return a `p x p` matrix.
#' @keywords internal
#' @export
model_info <- function(m, robust = FALSE) {
  if (!robust) return(m$OI)
  if (!is.null(m$robust_info)) return(m$robust_info)
  if (is.null(m$sandwich)) stop("This model has no robust covariance; use `robust = FALSE`")
  S <- solve(m$sandwich)
  (S + t(S)) / 2
}

#' Fit a model specification to a data set
#'
#' @param spec a specification from [mdpl_glm()], [mdpl_causl()] or
#'   [mdpl_custom()]; a formula is shorthand for `mdpl_glm(formula, family)`.
#' @param data a data frame.
#' @param family used only when `spec` is a formula.
#' @return an [mdpl_model][new_mdpl_model] object.
#' @export
mdpl_fit <- function(spec, data, family = stats::gaussian()) {
  spec <- as_spec(spec, family)
  m <- spec$fit(as.data.frame(data))
  if (!inherits(m, "mdpl_model")) stop("specification did not return an `mdpl_model`")
  m
}

as_spec <- function(spec, family = stats::gaussian()) {
  if (inherits(spec, "mdpl_spec")) return(spec)
  if (inherits(spec, "formula")) return(mdpl_glm(spec, family))
  stop("`spec` must be a specification from mdpl_glm()/mdpl_causl()/mdpl_custom() or a formula")
}

#' @export
print.mdpl_spec <- function(x, ...) {
  cat("<mdpl_spec> ", x$label, "\n", sep = "")
  invisible(x)
}

#' Parameters that determine the treatment effect
#'
#' @param fit an `mdpl_model`.
#' @param treat name of the treatment variable.
#' @return character vector of parameter names.
#' @export
treatment_terms <- function(fit, treat = "A") {
  tt <- fit$treatment_terms(treat)
  if (!length(tt)) stop("No parameter involves '", treat, "'")
  tt
}
