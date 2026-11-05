test_that("binomial fit matches glm and robust SEs match the sandwich", {
  set.seed(1)
  d <- sim_bin(500)
  f <- mdpl_fit(Y ~ A * X, d, binomial())
  g <- glm(Y ~ A * X, binomial(), d)
  expect_equal(f$par, coef(g), tolerance = 1e-8)
  expect_equal(unname(solve(f$OI)), unname(vcov(g)), tolerance = 1e-6)
  skip_if_not_installed("sandwich")
  expect_equal(unname(f$sandwich), unname(sandwich::sandwich(g)), tolerance = 1e-4)
})

test_that("log-likelihood matrix at the estimate equals glm's logLik", {
  set.seed(2)
  d <- sim_bin(300)
  f <- mdpl_fit(Y ~ A * X, d, binomial())
  g <- glm(Y ~ A * X, binomial(), d)
  expect_equal(sum(f$loglik_mat(matrix(f$par, 1))), as.numeric(logLik(g)), tolerance = 1e-8)
  fp <- mdpl_fit(Y ~ A * X, rpois_data <- transform(d, Y = rpois(nrow(d), 2)), poisson())
  gp <- glm(Y ~ A * X, poisson(), rpois_data)
  expect_equal(sum(fp$loglik_mat(matrix(fp$par, 1))), as.numeric(logLik(gp)), tolerance = 1e-8)
})

test_that("binary outcome: risk difference and log odds ratio", {
  set.seed(3)
  fit <- mdpl(Y ~ A * X, sim_bin(600), list(sim_bin(1500, 0), sim_bin(1500, 1)), family = binomial())
  nd <- data.frame(X = c(-1, 0, 1))
  rd <- mdpl_cate(fit, nd, trace = "t2", scale = "response")
  lo <- mdpl_cate(fit, nd, trace = "t2", scale = "link")
  th <- fit$posterior$t2$mean
  expect_equal(lo$estimate, unname(th["A"] + th["A:X"] * nd$X), tolerance = 1e-10)
  p1 <- plogis(th["(Intercept)"] + th["A"] + (th["X"] + th["A:X"]) * nd$X)
  p0 <- plogis(th["(Intercept)"] + th["X"] * nd$X)
  expect_equal(rd$estimate, unname(p1 - p0), tolerance = 1e-10)
  expect_true(all(rd$lower < rd$estimate & rd$estimate < rd$upper))
  ate <- mdpl_ate(fit, trace = "t2")
  expect_equal(nrow(ate), 1)
  b <- mdpl_cate_boot(fit, nd, B = 15, trace = "t2", seed = 1)
  expect_equal(b$estimate, rd$estimate, tolerance = 1e-10)
})

test_that("a hand-written logistic likelihood agrees with the GLM adapter", {
  set.seed(4)
  d <- sim_bin(400)
  ll <- function(theta, data) {
    eta <- theta[1] + theta[2] * data$A + theta[3] * data$X + theta[4] * data$A * data$X
    data$Y * eta - log1p(exp(eta))
  }
  contrast <- function(newdata, treat) cbind(0, 1, 0, newdata$X)
  spec <- mdpl_custom(ll, start = c(a = 0, A = 0, X = 0, AX = 0), contrast = contrast)
  m <- mdpl_fit(spec, d)
  g <- mdpl_fit(Y ~ A * X, d, binomial())
  expect_equal(unname(m$par), unname(g$par), tolerance = 1e-4)
  expect_equal(unname(m$OI), unname(g$OI), tolerance = 1e-3)
  expect_equal(unname(m$sandwich), unname(g$sandwich), tolerance = 1e-3)
  expect_equal(sum(m$loglik_mat(matrix(m$par, 1))), sum(g$loglik_mat(matrix(g$par, 1))), tolerance = 1e-6)

  fit <- mdpl(spec, d, list(sim_bin(800, 0)), treat = "A", target = c("A", "AX"))
  ct <- mdpl_cate(fit, data.frame(X = c(0, 1)), scale = "link")
  expect_equal(nrow(ct), 2)
})

test_that("mdpl_check_model flags a bad likelihood", {
  set.seed(5)
  d <- sim_bin(300)
  good <- mdpl_custom(function(theta, data) stats::dbinom(data$Y, 1, plogis(theta[1] + theta[2] * data$A), log = TRUE),
                      start = c(0, 0))
  expect_true(all(suppressMessages(utils::capture.output(r <- mdpl_check_model(good, d))) != "x"))
  expect_true(all(r$passed))
  bad <- mdpl_custom(function(theta, data) stats::dbinom(data$Y, 1, plogis(theta[1] + theta[2] * data$A + theta[3] * data$A), log = TRUE),
                     start = c(0, 0, 0))
  r2 <- suppressWarnings(utils::capture.output(rb <- mdpl_check_model(bad, d)))
  expect_false(all(rb$passed))
  expect_error(suppressWarnings(utils::capture.output(mdpl_check_model(bad, d, stop_on_fail = TRUE))))
})

test_that("a total-only likelihood works with taylor and cv(loss) but not elpd", {
  set.seed(5)
  ll_row <- function(theta, data) {
    eta <- theta[1] + theta[2] * data$A + theta[3] * data$X + theta[4] * data$A * data$X
    data$Y * eta - log1p(exp(eta))
  }
  ll_tot <- function(theta, data) sum(ll_row(theta, data))
  contrast <- function(newdata, treat) cbind(0, 1, 0, newdata$X)
  st <- c(a = 0, A = 0, X = 0, AX = 0)
  row <- mdpl_custom(ll_row, st, contrast = contrast)
  tot <- mdpl_custom(ll_tot, st, contrast = contrast, per_row = FALSE)
  d <- sim_bin(400)
  m <- mdpl_fit(tot, d)
  expect_null(m$loglik_mat)
  expect_null(m$sandwich)
  expect_equal(m$n, 400)
  expect_equal(unname(m$par), unname(mdpl_fit(row, d)$par), tolerance = 1e-4)

  obs <- list(sim_bin(800, 0))
  ft <- mdpl(tot, d, obs, treat = "A")
  expect_equal(ft$eta, mdpl(row, d, obs, treat = "A")$eta, tolerance = 1e-3)
  expect_error(mdpl(tot, d, obs, treat = "A", robust = TRUE), "robust")
  expect_error(mdpl(tot, d, obs, treat = "A", method = "elpd"), "per-row")
  expect_error(mdpl(tot, d, obs, treat = "A", method = "cv"), "loss")
  loss <- function(Theta, data) {
    eta <- Theta[, 1] + Theta[, 2] * data$A + Theta[, 3] * data$X + Theta[, 4] * data$A * data$X
    -(data$Y * eta - log1p(exp(eta)))
  }
  fc <- mdpl(tot, d, obs, treat = "A", method = "cv", loss = loss, folds = 4)
  expect_equal(dim(fc$eta), c(1, 1))
  expect_error(mdpl_custom(function(theta, data) rep(0, nrow(data)), st, per_row = FALSE) |>
                 mdpl_fit(d), "single total")
  chk <- capture.output(r <- mdpl_check_model(tot, d))
  expect_true("per-row checks skipped" %in% r$check)
})
