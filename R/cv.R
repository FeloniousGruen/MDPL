#' Select the borrowing weights by gold-standard cross-validation
#'
#' The reference criterion of estimator-level fusion. The trial is split into
#' folds; for each fold the trial model is refitted on the other folds and
#' fused with the observational sources, whose summaries do not change. The
#' fused estimate is scored on the held-out trial units. With the weights
#' \eqn{\eta}, fold \eqn{f} uses
#' \deqn{P^{(-f)} = I_e^{(-f)} + \gamma_f \sum_k \eta_k I_k, \quad
#'   \hat\theta^{(-f)} = (P^{(-f)})^{-1}\{I_e^{(-f)}\hat\theta_e^{(-f)} +
#'   \gamma_f \sum_k \eta_k I_k \hat\theta_k\},}
#' where \eqn{\gamma_f} is the fraction of the trial used for training. The
#' factor \eqn{\gamma_f} keeps the relative precision of the trial and the
#' sources as it will be in the full-data analysis.
#'
#' The validation loss defaults to the **trial loss**: minus the per-unit
#' log-likelihood (or partial-likelihood contribution) of the held-out units at
#' the fold's fused estimate. Supply `loss` to use another one; it is required
#' for models without per-row log-likelihoods.
#'
#' With `type = "if"` the refits are replaced by influence-function deletion,
#' \eqn{\hat\theta_e^{(-i)} \approx \hat\theta_e - \hat\psi_i / n}, which needs
#' the model to supply `dfbeta`, the `n x p` matrix of approximate changes in
#' the estimate when each unit is deleted ([gbi's][new_mdpl_model] survival
#' models do). Because the trial precision is rescaled by
#' \eqn{\gamma_i = (n-1)/n} together with the sources, the leave-one-out fused
#' estimate has the closed form
#' \eqn{\hat\theta_\eta - S_\eta^{-1} I_e\, \mathrm{dfbeta}_i}, so no model is
#' refitted.
#'
#' With `type = "press"` the trial model must be a linear regression on
#' fixed (pseudo-)outcomes: it supplies `linear = list(G, y, sigma2)`, the
#' design, the outcome and the residual variance (a Gaussian [mdpl_glm()] and
#' the `gbi` AIPW estimators do). The fused estimate is then the weighted
#' least-squares fit in which the trial rows have unit weight and the sources
#' enter with their information, and the exact leave-one-out squared
#' residuals follow from the PRESS identity \eqn{e_i/(1-h_i)}, so there is no
#' refitting and no influence-function approximation. The criterion is the
#' mean squared leave-one-out prediction error; `loss` is not used. It is
#' conditional on the outcomes, so for estimated pseudo-outcomes cross-fit
#' them.
#'
#' The weights are searched one source at a time over `grid`, starting from
#' zero (no borrowing), as in [elpd_eta()].
#'
#' @param fit_exp trial model.
#' @param fit_obs one model or a list of them.
#' @param spec,data_exp the specification and trial data, needed to refit the
#'   trial in each fold (`type = "kfold"`).
#' @param type `"kfold"` (exact refits), `"if"` (influence-function
#'   leave-one-out) or `"press"` (exact leave-one-out for a linear regression,
#'   see above).
#' @param folds number of folds, or an integer vector of fold labels (one per
#'   trial unit). Ignored for `type = "if"`.
#' @param loss optional validation loss, `function(Theta, data, model)`:
#'   `Theta` is a matrix with one row of parameters per row of `data` (the
#'   held-out trial units), `model` is the full trial model; it returns one
#'   loss per row, lower being better. `NULL` is the trial loss.
#' @param grid candidate weights for each source.
#' @param robust use the inverse sandwich covariance as the information.
#' @return an object of class `"cv_eta"` with `eta` (a `K x 1` matrix, column
#'   `cv`), `table` (candidate weight and loss for the last sweep of the first
#'   source, `tables` for all), `type`, `folds`.
#' @export
cv_eta <- function(fit_exp, fit_obs, spec = NULL, data_exp = NULL,
                   type = c("kfold", "if", "press"), folds = 10L, loss = NULL,
                   grid = c(seq(0, 0.2, 0.04), seq(0.25, 1, 0.05)), robust = FALSE) {
  type <- match.arg(type)
  if (inherits(fit_obs, "mdpl_model")) fit_obs <- list(fit_obs)
  if (type != "press" && is.null(loss) && is.null(fit_exp$loglik_mat))
    stop("this model has no per-row log-likelihood, so the default (trial-loss) validation is ",
         "unavailable; supply `loss`", call. = FALSE)
  K <- length(fit_obs); n <- fit_exp$n; nm <- names(fit_exp$par); p <- length(nm)
  Ik <- lapply(fit_obs, function(f) unname(model_info(f, robust)))
  thk <- lapply(fit_obs, function(f) unname(f$par))
  sum_obs <- function(eta, g = 1) {
    S <- matrix(0, p, p); b <- numeric(p)
    for (k in seq_len(K)) if (eta[k] > 0) {
      S <- S + g * eta[k] * Ik[[k]]
      b <- b + g * eta[k] * as.numeric(Ik[[k]] %*% thk[[k]])
    }
    list(S = S, b = b)
  }
  call_loss <- function(Theta, data) {
    dimnames(Theta) <- list(NULL, nm)
    out <- if ("model" %in% names(formals(loss)) || "..." %in% names(formals(loss)))
      loss(Theta, data, model = fit_exp) else loss(Theta, data)
    if (length(out) != nrow(Theta)) stop("`loss` must return one value per held-out row")
    out
  }

  if (type == "press") {
    lin <- fit_exp$linear
    if (is.null(lin)) stop("type = \"press\" needs a trial model that is a linear regression ",
                           "(it supplies `linear`): a Gaussian GLM or an AIPW estimator")
    if (!is.null(loss)) stop("`loss` is not used with type = \"press\"; it scores the squared ",
                             "leave-one-out prediction error")
    G <- as.matrix(lin$G); y <- as.numeric(lin$y)
    if (ncol(G) != p || nrow(G) != n) stop("`linear` does not match the trial model")
    GtG <- crossprod(G) / lin$sigma2; the <- unname(fit_exp$par)
    score <- function(eta) {
      o <- sum_obs(eta)
      Si <- solve(GtG + o$S)
      th <- Si %*% (GtG %*% the + o$b)
      e <- y - as.numeric(G %*% th)
      h <- rowSums((G %*% Si) * G) / lin$sigma2
      l <- (e / (1 - h))^2
      if (!all(is.finite(l))) Inf else mean(l)
    }
    fold_info <- NULL
  } else if (type == "if") {
    dfb <- fit_exp$dfbeta
    if (is.null(dfb) || nrow(dfb) != n)
      stop("type = \"if\" needs a model that supplies `dfbeta` (one row per trial unit)")
    Ie <- unname(model_info(fit_exp, robust)); the <- unname(fit_exp$par)
    if (!is.null(loss) && is.null(data_exp)) stop("A custom `loss` needs `data_exp`")
    dat <- if (!is.null(loss)) as.data.frame(data_exp)
    score <- function(eta) {
      o <- sum_obs(eta)
      S <- Ie + o$S
      th <- solve(S, as.numeric(Ie %*% the) + o$b)
      Th <- sweep(-dfb %*% t(solve(S, Ie)), 2, th, "+")
      l <- if (is.null(loss)) -diag(fit_exp$loglik_mat(Th)) else call_loss(Th, dat)
      if (!all(is.finite(l))) Inf else mean(l)
    }
    fold_info <- NULL
  } else {
    if (is.null(spec) || is.null(data_exp)) stop("type = \"kfold\" needs `spec` and `data_exp`")
    fid <- if (length(folds) == 1L) {
      if (folds < 2 || folds > n) stop("`folds` must be between 2 and the number of trial units")
      sample(rep_len(seq_len(folds), n))
    } else {
      if (length(folds) != n) stop("A vector `folds` needs one label per trial unit")
      as.integer(factor(folds))
    }
    data_exp <- as.data.frame(data_exp)
    subset_rows <- function(idx) {
      d <- data_exp[idx, , drop = FALSE]
      for (a in setdiff(names(attributes(data_exp)), names(attributes(d))))
        attr(d, a) <- attr(data_exp, a)
      d
    }
    fl <- lapply(sort(unique(fid)), function(f) {
      tr <- which(fid != f); te <- which(fid == f)
      m <- spec$fit(subset_rows(tr))
      if (!identical(names(m$par), nm))
        stop("The trial model refitted on a training fold has different parameters")
      list(te = te, g = length(tr) / n, Ie = unname(model_info(m, robust)),
           th = unname(m$par), test = data_exp[te, , drop = FALSE])
    })
    score <- function(eta) {
      tot <- 0
      for (f in fl) {
        o <- sum_obs(eta, f$g)
        th <- solve(f$Ie + o$S, as.numeric(f$Ie %*% f$th) + o$b)
        l <- if (is.null(loss)) -fit_exp$loglik_mat(matrix(th, 1L))[1L, f$te]
             else call_loss(matrix(th, length(f$te), p, byrow = TRUE), f$test)
        tot <- tot + sum(l)
      }
      if (is.finite(tot)) tot / n else Inf
    }
    fold_info <- fid
  }

  eta <- rep(0, K); tables <- vector("list", K)
  for (s in seq_len(if (K == 1L) 1L else 2L)) for (k in seq_len(K)) {
    v <- vapply(grid, function(g) { x <- eta; x[k] <- g; score(x) }, numeric(1))
    eta[k] <- grid[which.min(v)]
    tables[[k]] <- data.frame(eta = grid, loss = v)
  }
  eta <- matrix(eta, K, 1L, dimnames = list(names(fit_obs) %||% paste0("source", seq_len(K)), "cv"))
  structure(list(eta = eta, table = tables[[1]], tables = tables, K = K, type = type,
                 folds = fold_info), class = "cv_eta")
}

#' @export
print.cv_eta <- function(x, digits = 3, ...) {
  cat("Cross-validated borrowing weights (K =", x$K, "sources; ",
      switch(x$type, "if" = "influence-function leave-one-out", press = "PRESS leave-one-out", "k-fold"),
      ")\n", sep = "")
  print(round(x$eta, digits))
  invisible(x)
}
