reference_eta <- function(fit_exp, fit_obs, A, lam) {
  K <- length(fit_obs); phi_e <- fit_exp$par
  f <- function(x) {
    S <- fit_exp$OI
    for (k in 1:K) S <- S + x[k] * fit_obs[[k]]$OI
    Sinv <- solve(S)
    t1 <- 0
    for (k in 1:K) for (l in 1:K) {
      dk <- fit_obs[[k]]$par - phi_e; dl <- fit_obs[[l]]$par - phi_e
      t1 <- t1 + x[k] * x[l] *
        as.numeric(t(dk) %*% fit_obs[[k]]$OI %*% Sinv %*% A %*% Sinv %*% fit_obs[[l]]$OI %*% dl)
    }
    t1 + lam * sum(diag(A %*% Sinv))
  }
  optim(rep(0, K), f, method = "L-BFGS-B", lower = rep(0, K), upper = rep(1, K),
        control = list(maxit = 1000, factr = 1e-7, pgtol = 1e-8, fnscale = -1))
}

test_that("taylor_eta matches the double-sum reference at both traces", {
  set.seed(1)
  fe <- mdpl_fit(Y ~ A * X3, sim_data(400))
  fo <- lapply(c(0, 0.5, 1, 2), function(b) mdpl_fit(Y ~ A * X3, sim_data(800, b)))
  for (target in list(NULL, treatment_terms(fe))) {
    A <- mdpl_metric(fe, target)
    sel <- taylor_eta(fe, fo, A = A)
    for (lam in 1:2) {
      ref <- reference_eta(fe, fo, A, lam)
      expect_equal(unname(sel$eta[, lam]), ref$par, tolerance = 1e-3)
      expect_equal(sel$value[lam], ref$value, tolerance = 1e-6)
    }
  }
})

test_that("analytic gradient equals finite differences", {
  set.seed(2)
  fe <- mdpl_fit(Y ~ A * X3, sim_data(300))
  fo <- lapply(c(0.5, 1.5), function(b) mdpl_fit(Y ~ A * X3, sim_data(600, b)))
  a <- taylor_eta(fe, fo, target = treatment_terms(fe), gradient = TRUE)
  b <- taylor_eta(fe, fo, target = treatment_terms(fe), gradient = FALSE)
  expect_equal(a$eta, b$eta, tolerance = 1e-3)
})

test_that("two-trace weights are at least the one-trace weights", {
  set.seed(3)
  fe <- mdpl_fit(Y ~ A * X3, sim_data(300))
  fo <- lapply(c(0, 0, 0.3), function(b) mdpl_fit(Y ~ A * X3, sim_data(900, b)))
  sel <- taylor_eta(fe, fo)
  expect_true(all(sel$eta[, "t2"] >= sel$eta[, "t1"] - 1e-6))
})

test_that("Gaussian fit matches lm", {
  set.seed(4)
  d <- sim_data(300)
  f <- mdpl_fit(Y ~ A * X3, d)
  l <- lm(Y ~ A * X3, d)
  expect_equal(f$par, coef(l), tolerance = 1e-8)
  expect_equal(unname(f$OI), unname(solve(vcov(l)) * summary(l)$sigma^2 / f$phi), tolerance = 1e-6)
})

test_that("posterior: eta = 0 is the trial fit, eta = 1 adds the information", {
  set.seed(5)
  fe <- mdpl_fit(Y ~ A * X3, sim_data(300))
  fo <- mdpl_fit(Y ~ A * X3, sim_data(900))
  p0 <- mdpl_posterior(fe, fo, 0)
  expect_equal(p0$mean, fe$par)
  expect_equal(unname(p0$cov), unname(solve(fe$OI)), tolerance = 1e-10)
  p1 <- mdpl_posterior(fe, fo, 1)
  expect_equal(unname(p1$cov), unname(solve(fe$OI + fo$OI)), tolerance = 1e-10)
})

test_that("CATE and ATE are the treatment coefficients; several sources work", {
  set.seed(6)
  d <- sim_data(500)
  fit <- mdpl(Y ~ A * X3, d, list(a = sim_data(1000, 0), b = sim_data(1000, 1)))
  expect_equal(rownames(fit$eta), c("a", "b"))
  nd <- data.frame(X3 = c(0, 1, 2))
  ct <- mdpl_cate(fit, nd, trace = "t1")
  th <- fit$posterior$t1$mean
  expect_equal(ct$estimate, unname(th["A"] + th["A:X3"] * nd$X3), tolerance = 1e-10)
  expect_true(all(ct$lower < ct$estimate & ct$estimate < ct$upper))
  ate <- mdpl_ate(fit, trace = "t1")
  expect_equal(ate$estimate, mdpl_cate(fit, data.frame(X3 = mean(d$X3)), trace = "t1")$estimate)
  cf <- mdpl_coef(fit, trace = "t1")
  expect_equal(cf$term, c("A", "A:X3"))
})

test_that("pooled and elpd variants run", {
  set.seed(7)
  d <- sim_data(300)
  o <- list(sim_data(600, 0), sim_data(600, 0.5))
  fp <- mdpl(Y ~ A * X3, d, o, pooled = TRUE)
  expect_equal(nrow(fp$eta), 1)
  fe <- mdpl(Y ~ A * X3, d, o[[1]], method = "elpd", grid = c(0, 0.5, 1), n_draws = 200)
  expect_equal(colnames(fe$eta), "elpd")
  expect_true(fe$eta[1, 1] %in% c(0, 0.5, 1))
})

test_that("bootstrap summaries have the right shape and the point estimate is unchanged", {
  set.seed(8)
  fit <- mdpl(Y ~ A * X3, sim_data(200), sim_data(400, 0.5))
  nd <- data.frame(X3 = c(0, 1))
  b <- mdpl_cate_boot(fit, nd, B = 20, trace = "t2", seed = 1)
  expect_equal(nrow(b), 2)
  expect_equal(b$estimate, mdpl_cate(fit, nd, trace = "t2")$estimate)
  expect_equal(dim(attr(b, "draws")), c(20, 2, 1))
  ba <- mdpl_ate_boot(fit, B = 20, seed = 1)
  expect_equal(ba$trace, "t2")
  expect_equal(mdpl_ate_boot(fit, B = 20, seed = 1, trace = "t1")$trace, "t1")
  bc <- mdpl_coef_boot(fit, B = 20, trace = "t2", seed = 1)
  expect_equal(bc$term, c("A", "A:X3"))
})

test_that("bootstrap_CI = TRUE stores replicates and the effect functions use them", {
  set.seed(9)
  d <- sim_data(200); o <- sim_data(400, 0.5)
  fit0 <- mdpl(Y ~ A * X3, d, o)
  expect_null(fit0$boot)
  expect_error(mdpl_cate(fit0, data.frame(X3 = 1), ci = "bootstrap"))
  fit <- mdpl(Y ~ A * X3, d, o, bootstrap_CI = TRUE, B = 30, boot_seed = 1)
  expect_equal(dim(fit$boot$draws), c(30, 4, 2))
  nd <- data.frame(X3 = c(0, 1))
  auto <- mdpl_cate(fit, nd, trace = "t2")
  post <- mdpl_cate(fit, nd, trace = "t2", ci = "posterior")
  expect_equal(auto$estimate, post$estimate)
  expect_false(isTRUE(all.equal(auto$se, post$se)))
  expect_equal(auto$se, apply(fit$boot$draws[, , "t2"] %*% t(cbind(0, 1, 0, nd$X3)), 2, sd, na.rm = TRUE))
  expect_equal(nrow(mdpl_ate(fit, trace = "t2")), 1)
  expect_equal(mdpl_coef(fit, trace = "t2")$term, c("A", "A:X3"))
  b <- mdpl_cate_boot(fit, nd, trace = "t2")      # reuses the stored replicates
  expect_equal(b$se, auto$se)
})

test_that("t2 is the default trace everywhere and t1 stays available", {
  set.seed(10)
  d <- sim_data(200); o <- sim_data(400, 0.5)
  fit <- mdpl(Y ~ A * X3, d, o, bootstrap_CI = TRUE, B = 20, boot_seed = 1)
  nd <- data.frame(X3 = c(0, 1))
  expect_equal(mdpl_cate(fit, nd), mdpl_cate(fit, nd, trace = "t2"))
  expect_equal(mdpl_ate(fit), mdpl_ate(fit, trace = "t2"))
  expect_equal(mdpl_coef(fit), mdpl_coef(fit, trace = "t2"))
  expect_equal(mdpl_cate(fit, nd, ci = "posterior"), mdpl_cate(fit, nd, trace = "t2", ci = "posterior"))
  expect_equal(mdpl_cate_boot(fit, nd), mdpl_cate_boot(fit, nd, trace = "t2"))
  expect_equal(unique(mdpl_cate_boot(fit, nd)$trace), "t2")
  expect_false(isTRUE(all.equal(mdpl_cate(fit, nd, trace = "t1")$estimate,
                                mdpl_cate(fit, nd)$estimate)))
})

test_that("target is optional, defaults to all parameters, and only applies to taylor", {
  set.seed(11)
  d <- sim_data(200); o <- sim_data(400, 0.5)
  fit <- mdpl(Y ~ A * X3, d, o)
  expect_null(fit$target)
  expect_equal(fit$eta, mdpl(Y ~ A * X3, d, o, target = "all")$eta)
  fc <- mdpl(Y ~ A * X3, d, o, target = "cate")
  expect_equal(fc$target, c("A", "A:X3"))
  expect_false(isTRUE(all.equal(fc$eta, fit$eta)))
  expect_warning(mdpl(Y ~ A * X3, d, o, method = "elpd", target = "cate"), "only applies")
  expect_no_warning(mdpl(Y ~ A * X3, d, o, method = "elpd"))
})

test_that("sources that target different parameters cannot be merged", {
  set.seed(6)
  fe <- mdpl_fit(Y ~ A * X3, sim_data(300))
  fo <- mdpl_fit(Y ~ A + X3, sim_data(300))
  expect_error(mdpl_posterior(fe, fo, 0.5), "do not|cannot be merged")
  expect_true(check_same_target(fe, mdpl_fit(Y ~ A * X3, sim_data(300))))
})

test_that("mdpl() returns the fit, weights and coef / ATE / CATE estimates", {
  set.seed(7)
  nd <- data.frame(X3 = c(0, 1))
  f <- mdpl(Y ~ A * X3, sim_data(300), list(sim_data(600)), cate = nd)
  expect_named(f$estimates, c("trace", "coef", "ate", "cate", "cate_at"))
  expect_equal(f$estimates$coef, mdpl_coef(f))
  expect_equal(f$estimates$ate, mdpl_ate(f))
  expect_equal(f$estimates$cate_at, mdpl_cate(f, nd))
  expect_s3_class(f$estimates$cate, "mdpl_cate_fn")
  expect_equal(f$estimates$cate(nd, scale = "link"), mdpl_cate(f, nd))
  cf <- attr(f$estimates$cate, "coefficients")
  expect_equal(cf$covariate, c("1", "X3"))
  expect_match(attr(f$estimates$cate, "text"), "^CATE = .*\\*X3")
  expect_output(print(f), "CATE = ")
  g <- mdpl(Y ~ A * X3, sim_data(300), list(sim_data(600)), coef = FALSE, cate = FALSE)
  expect_named(g$estimates, c("trace", "ate"))
  h <- mdpl(Y ~ A * X3, sim_data(300), list(sim_data(600)), cate = FALSE, ate = FALSE, coef = FALSE)
  expect_null(h$estimates)
})

test_that("GLM offsets are used", {
  set.seed(20)
  n <- 3000; A <- rbinom(n, 1, .5); fu <- runif(n, 0.5, 4)
  d <- data.frame(A = A, fu = fu, y = rpois(n, fu * exp(-1 + 0.4 * A)))
  m <- mdpl_fit(mdpl_glm(y ~ A + offset(log(fu)), poisson()), d)
  expect_equal(unname(m$par), c(-1, 0.4), tolerance = 0.1)
  ref <- glm(y ~ A + offset(log(fu)), poisson(), d)
  expect_equal(unname(m$par), unname(coef(ref)), tolerance = 1e-6)
  fit <- mdpl(y ~ A + offset(log(fu)), d, list(d), treat = "A", family = poisson())
  expect_equal(mdpl_coef(fit)$estimate, 0.4, tolerance = 0.1)
})
