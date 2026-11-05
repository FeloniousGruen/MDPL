# Coverage of posterior vs bootstrap intervals for the CATE at X = 1.
# One trial (n = 300), one observational source (n = 1000) with unmeasured
# confounding of strength `bias`; 300 replicates, B = 200 bootstrap replicates
# each, about 13 minutes on 12 cores. Truth is 0.2 + 0.1 * 1 = 0.3.
#   Rscript coverage-study.R
library(MDPL)
library(parallel)

sim <- function(n, bias = 0) {
  X <- rnorm(n, 1); U <- rnorm(n)
  A <- rbinom(n, 1, if (bias == 0) 0.5 else plogis(X - 1 + bias * U))
  data.frame(X = X, A = A, Y = 0.5 + X + 0.2 * A + 0.1 * A * X + bias * U + rnorm(n))
}
nd <- data.frame(X = 1); truth <- 0.3
inside <- function(z) z$lower <= truth & truth <= z$upper

one <- function(i, bias) {
  set.seed(1000 * bias + i)
  f <- mdpl(Y ~ A * X, sim(300), list(sim(1000, bias)), treat = "A",
            bootstrap_CI = TRUE, B = 200, boot_seed = i)
  do.call(rbind, lapply(c("t1", "t2"), function(tr) {
    p <- mdpl_cate(f, nd, trace = tr, ci = "posterior")
    b <- mdpl_cate(f, nd, trace = tr, ci = "bootstrap")
    data.frame(bias = bias, trace = tr, seed = i,
               post_cover = inside(p), boot_cover = inside(b),
               post_width = p$upper - p$lower, boot_width = b$upper - b$lower)
  }))
}

res <- do.call(rbind, lapply(c(0, 0.5, 1, 2), function(b)
  do.call(rbind, mclapply(1:300, one, bias = b, mc.cores = 12))))
aggregate(cbind(post_cover, boot_cover, post_width, boot_width) ~ bias + trace, res, mean)
