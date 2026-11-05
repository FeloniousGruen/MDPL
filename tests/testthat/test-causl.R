test_that("causl frugal adapter: CATE equals the coefficient contrast", {
  skip_if_not_installed("causl")
  set.seed(3)
  forms <- list(Z ~ C1, Y ~ X * C1, ~ 1)
  fam <- list(1, c(5, 1), 1, 1)
  pars <- list(C1 = list(beta = 0, phi = 1), Z = list(beta = c(0, 0.5), phi = 1),
               X = list(beta = 0), Y = list(beta = c(0.5, 0.8, 0.2, 0.3), phi = 1),
               cop = list(beta = matrix(0.5, 1, 1)))
  f <- list(Z ~ C1, list(X ~ 1, C1 ~ 1), Y ~ X * C1, ~ 1)
  d_e <- causl::rfrugalParam(n = 200, formulas = f, family = fam, pars = pars)
  d_o <- causl::rfrugalParam(n = 400, formulas = f, family = fam, pars = pars)
  spec <- mdpl_causl(forms)
  fe <- mdpl_fit(spec, d_e)
  expect_true(all(c("Y:X", "Y:C1:X") %in% names(fe$par)))
  expect_equal(treatment_terms(fe, "X"), c("Y:X", "Y:C1:X"))

  fit <- mdpl(spec, d_e, list(d_o), treat = "X", method = "elpd",
              grid = c(0, 0.5, 1), n_draws = 100)
  nd <- data.frame(C1 = c(-1, 0, 2))
  ct <- mdpl_cate(fit, nd)
  th <- fit$posterior$elpd$mean
  expect_equal(ct$estimate, unname(th["Y:X"] + th["Y:C1:X"] * nd$C1), tolerance = 1e-8)
  expect_equal(mdpl_ate(fit)$estimate, mean(th["Y:X"] + th["Y:C1:X"] * d_e$C1), tolerance = 1e-8)

  ft <- mdpl(spec, d_e, list(d_o), treat = "X")
  expect_equal(colnames(ft$eta), c("t1", "t2"))
})
