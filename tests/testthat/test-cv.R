test_that("cv_eta down-weights a biased source and keeps a clean one", {
  set.seed(11)
  d <- sim_data(400)
  spec <- mdpl_glm(Y ~ A * X3)
  fe <- mdpl_fit(spec, d)
  fo <- list(clean = mdpl_fit(spec, sim_data(1500)), biased = mdpl_fit(spec, sim_data(1500, 2)))
  sel <- cv_eta(fe, fo, spec = spec, data_exp = d, folds = 5, robust = TRUE)
  expect_s3_class(sel, "cv_eta")
  expect_equal(dim(sel$eta), c(2L, 1L))
  expect_true(all(sel$eta >= 0 & sel$eta <= 1))
  expect_gt(sel$eta["clean", "cv"], sel$eta["biased", "cv"])
})

test_that("a custom loss is used and gets one parameter row per held-out unit", {
  set.seed(12)
  d <- sim_data(200)
  spec <- mdpl_glm(Y ~ A * X3)
  fe <- mdpl_fit(spec, d)
  fo <- list(mdpl_fit(spec, sim_data(600)))
  seen <- 0
  l <- function(Theta, data, model) {
    seen <<- seen + 1
    expect_equal(nrow(Theta), nrow(data))
    expect_true(inherits(model, "mdpl_model"))
    (data$Y - Theta[, "(Intercept)"] - Theta[, "X3"] * data$X3 - Theta[, "A"] * data$A -
       Theta[, "A:X3"] * data$A * data$X3)^2
  }
  sel <- cv_eta(fe, fo, spec = spec, data_exp = d, folds = 4, loss = l, grid = c(0, 0.5, 1))
  expect_gt(seen, 0)
  expect_true(sel$eta[1, 1] %in% c(0, 0.5, 1))
  expect_error(cv_eta(fe, fo, spec = spec, data_exp = d, loss = function(Theta, data) 1),
               "one value per")
})

test_that("method = 'cv' and a user metric run through mdpl()", {
  set.seed(13)
  d <- sim_data(250)
  fit <- mdpl(Y ~ A * X3, d, list(sim_data(800)), method = "cv", folds = 5)
  expect_equal(colnames(fit$eta), "cv")
  expect_match(capture.output(print(fit))[2], "held-out trial units")
  H <- diag(length(fit$fit_exp$par))
  fit2 <- mdpl(Y ~ A * X3, d, list(sim_data(800)), metric = H)
  expect_match(capture.output(print(fit2))[2], "user-supplied metric")
  expect_error(mdpl(Y ~ A * X3, d, list(sim_data(800)), metric = diag(2)), "p x p")
  expect_error(cv_eta(fit$fit_exp, fit$fit_obs, type = "if"), "dfbeta")
})

test_that("type = 'press' equals brute-force leave-one-out of the joint weighted fit", {
  set.seed(21)
  d <- sim_data(80)
  spec <- mdpl_glm(Y ~ A * X3)
  fe <- mdpl_fit(spec, d)
  fo <- list(mdpl_fit(spec, sim_data(500)))
  grid <- c(0, 0.3, 1)
  sel <- cv_eta(fe, fo, type = "press", grid = grid)
  lin <- fe$linear
  G <- lin$G; y <- lin$y; s2 <- lin$sigma2
  Ik <- unname(fo[[1]]$OI); thk <- unname(fo[[1]]$par)
  brute <- vapply(grid, function(g) {
    mean(vapply(seq_len(nrow(G)), function(i) {
      Gi <- G[-i, , drop = FALSE]
      th <- solve(crossprod(Gi) / s2 + g * Ik, crossprod(Gi, y[-i]) / s2 + g * Ik %*% thk)
      (y[i] - sum(G[i, ] * th))^2
    }, 1))
  }, 1)
  expect_equal(sel$table$loss, brute)
  expect_error(cv_eta(mdpl_fit(mdpl_glm(Y ~ A, binomial()), transform(d, Y = as.numeric(Y > 0))),
                      fo, type = "press"), "linear regression")
})
