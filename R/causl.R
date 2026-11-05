#' Specification: causl frugal (copula) model
#'
#' Wraps the frugal parametrisation of the 'causl' package. The joint law is
#' split into four blocks: the pre-treatment variables \eqn{Z}, the treatment
#' \eqn{X}, the outcome margin \eqn{q(y \mid do(x), c)}, and a copula tying the
#' margins together. Because the treatment model does not share parameters
#' with the others, it is left out of the fit; only the \eqn{Z} block, the
#' outcome block and the copula are estimated, and the treatment effect is the
#' contrast of the outcome margin's linear predictor at treatment 1 versus 0
#' (identity link).
#'
#' @param formulas a flat list of formulas in fitting order: one per
#'   pre-treatment variable (the \eqn{Z} block), then the outcome formula,
#'   then the copula formula (`~ ...`). The treatment appears only on the
#'   right-hand side of the outcome (and possibly the copula) formula.
#' @param family numeric family codes, one per formula including the copula
#'   (last entry). Supported margins: `1` (Gaussian) and `5` (binary). The
#'   copula family is the last entry (`1` = Gaussian).
#' @param outcome name of the outcome (left-hand side of its formula).
#' @param kwd left-hand side keyword of the copula formula.
#'
#' @details Needs the 'causl' package (a suggested dependency). With causl
#' 0.10.4.9000 a copula formula other than `~ 1` (for example `~ X * C1`)
#' fails inside causl when there is more than one copula pair, that is, with
#' several `Z` variables; use `~ 1` for the copula in that case.
#' @return an object of class `"mdpl_spec"`.
#' @examples
#' \dontrun{
#' spec <- mdpl_causl(list(Z ~ C1, Y ~ X * C1, ~ X * C1))
#' }
#' @export
mdpl_causl <- function(formulas, family = rep(1, length(formulas)), outcome = "Y", kwd = "cop") {
  if (!requireNamespace("causl", quietly = TRUE))
    stop("Package 'causl' is needed for mdpl_causl(); install it first")
  info <- causl_info(formulas, family, outcome, kwd)
  structure(list(fit = function(data) fit_causl_model(info, as.data.frame(data)),
                 label = paste0("causl frugal model, margins: ", paste(info$vars, collapse = ", ")),
                 info = info), class = "mdpl_spec")
}

causl_info <- function(formulas, family, outcome, kwd) {
  nf <- length(formulas)
  if (length(family) != nf) stop("`family` needs one code per formula, including the copula")
  forms <- causl:::tidy_formulas(formulas, kwd = kwd)
  full_form <- causl:::merge_formulas(forms)
  vars <- causl:::lhs(forms[-nf])
  if (!outcome %in% vars) stop("`outcome` must be the left-hand side of one of the formulas")
  nv <- nf - 1L
  ncop <- choose(nv, 2)
  beta_m <- matrix(0, nrow = max(unlist(full_form$wh)), ncol = nv + ncop)
  for (i in seq_len(nv)) beta_m[full_form$wh[[i]], i] <- 1
  for (j in seq_len(ncop)) beta_m[full_form$wh[[nf]], nv + j] <- 1
  fam <- lapply(family[seq_len(nv)], function(f) {
    if (f == 1) causl::gaussian_causl_fam()
    else if (f == 5) causl::binomial_causl_fam()
    else stop("only Gaussian (1) and binary (5) margins are supported")
  })
  list(formulas = formulas, codes = family, forms = forms, full_form = full_form,
       masks = list(beta_m = beta_m, phi_m = as.numeric(family[seq_len(nv)] == 1)),
       vars = vars, family = fam, fam_cop = family[[nf]], inCop = seq_len(nv),
       outcome = outcome)
}

causl_par_names <- function(info, mm) {
  bm <- info$masks$beta_m
  nv <- length(info$vars)
  cols <- c(info$vars, paste0("cop", seq_len(ncol(bm) - nv)))
  nm <- character(0)
  for (j in seq_len(ncol(bm)))
    nm <- c(nm, paste0(cols[j], ":", colnames(mm)[which(bm[, j] > 0)]))
  c(nm, paste0("phi:", info$vars[info$masks$phi_m > 0]))
}

fit_causl_model <- function(info, data) {
  fit <- causl::fit_causl(data, formulas = info$formulas, family = info$codes)
  nm <- causl_par_names(info, fit$mm)
  par <- stats::setNames(fit$par, nm)
  OI <- fit$OI
  sand <- fit$sandwich
  dat <- data[, info$vars, drop = FALSE]
  mm <- fit$mm
  msks <- info$masks
  np <- sum(msks$beta_m > 0)
  loglik_mat <- function(Theta) {
    Theta <- matrix(Theta, ncol = length(par))
    t(vapply(seq_len(nrow(Theta)), function(s) {
      m2 <- msks
      m2$beta_m[msks$beta_m > 0] <- Theta[s, seq_len(np)]
      m2$phi_m[msks$phi_m > 0] <- Theta[s, -seq_len(np)]
      as.vector(causl:::ll(dat, mm, m2$beta_m, phi = m2$phi_m, inCop = info$inCop,
                           family = info$family, fam_cop = info$fam_cop))
    }, numeric(nrow(dat))))
  }
  contrast <- function(newdata, treat, scale = "response") {
    nd1 <- nd0 <- as.data.frame(newdata)
    nd1[[treat]] <- 1; nd0[[treat]] <- 0
    D <- stats::model.matrix(info$full_form$formula, nd1) -
      stats::model.matrix(info$full_form$formula, nd0)
    bm <- msks$beta_m
    y <- match(info$outcome, info$vars)
    rows <- which(bm[, y] > 0)
    offset <- sum(bm[, seq_len(y - 1), drop = FALSE] > 0)
    L <- matrix(0, nrow(D), length(par))
    L[, offset + seq_along(rows)] <- D[, rows, drop = FALSE]
    L
  }
  effect <- function(Theta, newdata, treat, scale = "response")
    matrix(Theta, ncol = length(par)) %*% t(contrast(newdata, treat))
  new_mdpl_model(par, OI, sand, nrow(data), loglik_mat, effect, contrast,
                 treatment_terms = function(treat) {
                   pieces <- strsplit(names(par), ":", fixed = TRUE)
                   hit <- vapply(pieces, function(z) z[1] == info$outcome && treat %in% z[-1], logical(1))
                   names(par)[hit]
                 },
                 fit = fit, info = info, label = "causl")
}
