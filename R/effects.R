# t2 (leave-one-out) is the default trace; fall back to the last one if absent
default_trace <- function(x) {
  nm <- names(x$posterior)
  if ("t2" %in% nm) "t2" else nm[length(nm)]
}

as_posterior <- function(x, trace = NULL) {
  if (inherits(x, "mdpl_posterior")) return(x)
  if (inherits(x, "mdpl")) {
    trace <- trace %||% default_trace(x)
    if (!trace %in% names(x$posterior))
      stop("`trace` must be one of: ", paste(names(x$posterior), collapse = ", "))
    return(x$posterior[[trace]])
  }
  stop("Expected an `mdpl` or `mdpl_posterior` object")
}

prep_newdata <- function(newdata, treat) {
  newdata <- as.data.frame(newdata)
  if (!treat %in% names(newdata)) newdata[[treat]] <- 0
  newdata
}

# S x m matrix of effects for the S parameter vectors in Theta (m = 1 for the ATE)
effect_matrix <- function(model, Theta, newdata, treat, scale, ate) {
  Theta <- matrix(Theta, ncol = length(model$par))
  L <- model$contrast(newdata, treat, scale)
  E <- if (!is.null(L)) Theta %*% t(L) else model$effect(Theta, newdata, treat, scale)
  if (ate) E <- matrix(rowMeans(E), ncol = 1)
  E
}

effect_summary <- function(post, newdata, treat, scale, level, ate, n_draws) {
  m <- post$model
  nd <- prep_newdata(newdata, treat)
  L <- m$contrast(nd, treat, scale)
  a <- (1 - level) / 2
  if (!is.null(L)) {
    if (ate) L <- matrix(colMeans(L), nrow = 1)
    est <- as.numeric(L %*% post$mean)
    se <- sqrt(rowSums((L %*% post$cov) * L))
    z <- stats::qnorm(1 - a)
    return(data.frame(estimate = est, se = se, lower = est - z * se, upper = est + z * se))
  }
  est <- effect_matrix(m, post$mean, nd, treat, scale, ate)[1, ]
  dr <- effect_matrix(m, mdpl_draws(post, n_draws), nd, treat, scale, ate)
  data.frame(estimate = est, se = apply(dr, 2, stats::sd),
             lower = apply(dr, 2, stats::quantile, a, names = FALSE),
             upper = apply(dr, 2, stats::quantile, 1 - a, names = FALSE))
}

resolve_ci <- function(x, ci) {
  has <- inherits(x, "mdpl") && !is.null(x$boot)
  if (ci == "auto") return(if (has) "bootstrap" else "posterior")
  if (ci == "bootstrap" && !has)
    stop("No stored bootstrap: refit with `mdpl(..., bootstrap_CI = TRUE)` or use mdpl_cate_boot()")
  ci
}

# bootstrap summary from the replicates stored in an `mdpl` fit
stored_boot_summary <- function(x, trace, fun, level, interval) {
  trace <- trace %||% default_trace(x)
  bt <- list(draws = x$boot$draws[, , trace, drop = FALSE], traces = trace)
  finish_boot(x, bt, fun, level, interval, NULL)[, c("estimate", "se", "lower", "upper")]
}

resolve_treat <- function(post, treat) {
  treat <- treat %||% post$treat
  if (is.null(treat)) stop("Please give `treat`")
  treat
}

#' Conditional average treatment effect
#'
#' The effect of switching the treatment from 0 to 1 at each row of `newdata`,
#' with its posterior standard error and credible interval.
#'
#' For models where the effect is linear in the parameters (identity link, or
#' `scale = "link"`) the summary is exact. Otherwise (e.g. the risk difference
#' of a logistic model) the effect is computed for `n_draws` posterior draws:
#' the estimate is the plug-in value at the posterior mean and the interval
#' uses the draw quantiles.
#'
#' @param x an `mdpl` object or an `mdpl_posterior`.
#' @param newdata data frame of covariate values (the treatment column is
#'   overwritten).
#' @param treat name of the treatment; defaults to the one used in the fit.
#' @param scale `"response"` (difference of the mean outcome; risk difference
#'   for a binary outcome) or `"link"` (difference of the linear predictor,
#'   e.g. a log odds ratio).
#' @param level credible level.
#' @param trace which trace multiplier to use for an `mdpl` object. The default
#'   is `"t2"` (leave-one-out); `"t1"` is the more conservative Jensen version
#'   and is available on request.
#' @param n_draws posterior draws used for nonlinear effects.
#' @param ci `"auto"` (default) uses the stored bootstrap if the fit was made
#'   with `bootstrap_CI = TRUE`, otherwise the posterior; or force
#'   `"posterior"` / `"bootstrap"`.
#' @param interval for bootstrap intervals, `"percentile"` or `"normal"`.
#' @return `newdata` with columns `estimate`, `se`, `lower`, `upper`.
#' @export
mdpl_cate <- function(x, newdata, treat = NULL, scale = c("response", "link"),
                      level = 0.95, trace = NULL, n_draws = 4000L,
                      ci = c("auto", "posterior", "bootstrap"),
                      interval = c("percentile", "normal")) {
  scale <- match.arg(scale); interval <- match.arg(interval)
  ci <- resolve_ci(x, match.arg(ci))
  post <- as_posterior(x, trace)
  treat <- resolve_treat(post, treat)
  nd <- as.data.frame(newdata)
  res <- if (ci == "bootstrap") {
    ndt <- prep_newdata(nd, treat)
    stored_boot_summary(x, trace, function(Th) effect_matrix(x$fit_exp, Th, ndt, treat, scale, FALSE),
                        level, interval)
  } else effect_summary(post, nd, treat, scale, level, ate = FALSE, n_draws)
  cbind(nd[, setdiff(names(nd), treat), drop = FALSE], res)
}

#' Average treatment effect
#'
#' The average of the CATE over the covariate distribution in `data` (the
#' trial by default), i.e. the standardised effect.
#'
#' @inheritParams mdpl_cate
#' @param data data frame whose covariate distribution defines the target
#'   population; default the trial data of an `mdpl` fit.
#' @return a one-row data frame with `estimate`, `se`, `lower`, `upper`.
#' @export
mdpl_ate <- function(x, data = NULL, treat = NULL, scale = c("response", "link"),
                     level = 0.95, trace = NULL, n_draws = 4000L,
                     ci = c("auto", "posterior", "bootstrap"),
                     interval = c("percentile", "normal")) {
  scale <- match.arg(scale); interval <- match.arg(interval)
  ci <- resolve_ci(x, match.arg(ci))
  post <- as_posterior(x, trace)
  treat <- resolve_treat(post, treat)
  if (is.null(data)) data <- if (inherits(x, "mdpl")) x$data_exp else stop("Please give `data`")
  if (ci == "bootstrap") {
    ndt <- prep_newdata(data, treat)
    return(stored_boot_summary(x, trace, function(Th) effect_matrix(x$fit_exp, Th, ndt, treat, scale, TRUE),
                               level, interval))
  }
  effect_summary(post, data, treat, scale, level, ate = TRUE, n_draws)
}

# Outcome-block parameters of a causl model are named "<outcome>:<term>"; show
# them as plain regression terms.
display_terms <- function(model, terms) {
  block <- model$info$outcome
  if (is.null(block)) return(terms)
  pre <- paste0(block, ":")
  ifelse(startsWith(terms, pre), substring(terms, nchar(pre) + 1L), terms)
}

# Accept either internal parameter names or the displayed (prefix-free) ones.
resolve_terms <- function(model, terms, par_names) {
  block <- model$info$outcome
  if (is.null(block)) return(terms)
  alt <- paste0(block, ":", terms)
  ifelse(terms %in% par_names, terms, ifelse(alt %in% par_names, alt, terms))
}

#' Posterior summary of the treatment-effect parameters
#'
#' @param x an `mdpl` object or `mdpl_posterior`.
#' @param terms parameter names; default the parameters that determine the
#'   treatment effect.
#' @param level credible level.
#' @param trace,ci,interval as in [mdpl_cate()].
#' @return a data frame with `term`, `estimate`, `se`, `lower`, `upper`.
#' @export
mdpl_coef <- function(x, terms = NULL, level = 0.95, trace = NULL,
                      ci = c("auto", "posterior", "bootstrap"),
                      interval = c("percentile", "normal")) {
  interval <- match.arg(interval)
  ci <- resolve_ci(x, match.arg(ci))
  post <- as_posterior(x, trace)
  terms <- terms %||% treatment_terms(post$model, post$treat %||% "A")
  terms <- resolve_terms(post$model, terms, names(post$mean))
  shown <- display_terms(post$model, terms)
  if (ci == "bootstrap") {
    ix <- match(terms, names(post$mean))
    r <- stored_boot_summary(x, trace, function(Th) matrix(Th, ncol = length(post$mean))[, ix, drop = FALSE],
                             level, interval)
    return(cbind(term = shown, r))
  }
  est <- post$mean[terms]
  se <- sqrt(diag(post$cov)[terms])
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(term = shown, estimate = unname(est), se = unname(se),
             lower = unname(est - z * se), upper = unname(est + z * se))
}

#' The CATE as a function (and formula) of the covariates
#'
#' For a model in which the treatment effect is linear in the parameters on the
#' link scale, the CATE is a formula in the covariates: the treatment
#' coefficient plus each treatment-by-covariate coefficient times that
#' covariate (for example `0.20 + 0.10*X` for `Y ~ A * X`). This returns that
#' formula together with a function that evaluates it, with intervals, on any
#' data frame (bootstrap intervals if the fit stored a bootstrap, otherwise
#' posterior intervals).
#'
#' The formula is on the link scale, which is the response scale for an
#' identity link (a mean difference) and a log odds or log rate ratio for
#' logistic and Poisson models. The returned function can also evaluate on the
#' response scale (for example a risk difference) through its `scale` argument,
#' but then there is no closed formula. Models for which the effect is not
#' linear in the parameters (some [mdpl_custom()] and [mdpl_causl()] models)
#' have no formula; use [mdpl_cate()] there.
#'
#' @param x an `mdpl` object (the formula needs the trial data, so an
#'   `mdpl_posterior` alone is not enough).
#' @param treat,trace,level as in [mdpl_cate()].
#' @return an object of class `"mdpl_cate_fn"`, a function
#'   `f(newdata, scale = "link", level = 0.95)` returning `newdata` with
#'   `estimate`, `se`, `lower`, `upper`. `attr(f, "coefficients")`,
#'   `attr(f, "text")` and `attr(f, "scale")` hold the table of terms, the
#'   formula as text and the scale.
#' @examples
#' set.seed(1)
#' d <- function(n) { X <- rnorm(n); A <- rbinom(n, 1, .5)
#'   data.frame(X = X, A = A, Y = 0.2 * A + 0.1 * A * X + X + rnorm(n)) }
#' fit <- mdpl(Y ~ A * X, d(300), list(d(800)))
#' f <- mdpl_cate_formula(fit)
#' f
#' f(data.frame(X = c(-1, 0, 1)))
#' @export
mdpl_cate_formula <- function(x, treat = NULL, trace = NULL, level = 0.95) {
  if (!inherits(x, "mdpl")) stop("`x` must be an `mdpl` fit")
  post <- as_posterior(x, trace)
  treat <- resolve_treat(post, treat)
  m <- post$model
  L <- tryCatch(m$contrast(prep_newdata(x$data_exp, treat), treat, "link"),
                error = function(e) NULL)
  if (is.null(L)) stop("the treatment effect of this model is not linear in its parameters, ",
                       "so it has no formula; use mdpl_cate() at covariate values")
  tt <- treatment_terms(m, treat)
  extra <- setdiff(colnames(L)[colSums(abs(L)) > 0], tt)
  if (length(extra))
    stop("the treatment effect depends on parameters that are not treatment terms: ",
         paste(extra, collapse = ", "))
  est <- post$mean[tt]; se <- sqrt(diag(post$cov)[tt])
  z <- stats::qnorm(1 - (1 - level) / 2)
  block <- m$info$outcome
  cov_part <- vapply(strsplit(tt, ":", fixed = TRUE), function(s) {
    if (!is.null(block) && s[1] == block) s <- s[-1]
    s <- setdiff(s, treat); if (length(s)) paste(s, collapse = ":") else "1"
  }, "")
  tab <- data.frame(term = display_terms(m, tt), covariate = cov_part, estimate = unname(est), se = unname(se),
                    lower = unname(est - z * se), upper = unname(est + z * se))
  if (resolve_ci(x, "auto") == "bootstrap") {
    b <- mdpl_coef(x, tt, level = level, trace = trace, ci = "bootstrap")
    tab[, c("se", "lower", "upper")] <- b[, c("se", "lower", "upper")]
  }
  fmt <- function(v) formatC(v, digits = 3, format = "g")
  pieces <- ifelse(cov_part == "1", fmt(tab$estimate),
                   paste0(fmt(abs(tab$estimate)), "*", cov_part))
  txt <- pieces[1]
  for (i in seq_along(pieces)[-1])
    txt <- paste0(txt, if (tab$estimate[i] < 0 && tab$covariate[i] != "1") " - " else " + ",
                  pieces[i])
  txt <- paste0("CATE = ", txt)
  tr <- trace %||% default_trace(x)
  fn <- function(newdata, scale = c("link", "response"), level = 0.95, n_draws = 4000L) {
    mdpl_cate(x, newdata, treat = treat, scale = match.arg(scale), level = level,
              trace = tr, n_draws = n_draws)
  }
  structure(fn, class = c("mdpl_cate_fn", "function"), coefficients = tab, text = txt,
            scale = "link", treat = treat, trace = tr)
}

#' @export
print.mdpl_cate_fn <- function(x, digits = 3, ...) {
  cat(attr(x, "text"), "   [link scale, ", attr(x, "trace"), "]\n", sep = "")
  tab <- attr(x, "coefficients")
  print(format(tab[, c("term", "estimate", "se", "lower", "upper")], digits = digits),
        row.names = FALSE)
  cat("Call it on a data frame of covariates to get the CATE with intervals.\n")
  invisible(x)
}
