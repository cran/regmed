## Tests for regmed


context("Testing regmed")

expect2 <- c("med.1", "med.2")
expect5 <- c("med.1", "med.2", "med.74", "med.88", "med.99")
load("regmedfits.RData")

data(medsim)
filter2 <- regmed.prefilter(x[,1], med, y[,1], k=2)

filter5 <- regmed.prefilter(x[,1], med, y[,1], k=5)
test_that("prefilter", {
  expect_equal(colnames(filter2$med), expected=expect2)
  expect_equal(colnames(filter5$med), expected=expect5)
  })


fit.grid5 <- regmed.grid(filter5$x, filter5$med, filter5$y, lambda.vec=seq(1, .1, by=-.1), frac.lasso=.8)
test_that("regmed_grid", {
  expect_equal(fit.grid5$grid.data$df, grid5$grid.data$df)
  expect_equal(fit.grid5$grid.data$bic, expected=grid5$grid.data$bic, tol=1e-2)
  })

fit.lam4 <- regmed.fit(filter5$x, filter5$med, filter5$y, lambda=.4, frac.lasso=.8)

edges4 <- regmed.edges(fit.lam4)
edgescoeff <- c(0.3579, 0.08275, 0.5104) # edges4$edges$coeff

#lam4 <- fit.lam4; grid5 <- fit.grid5
#save(lam4, grid5, file="regmedfits.RData")     

test_that("regmed_single", {
  expect_equal(fit.lam4$alpha[,1], lam4$alpha[,1], tol=1e-4)
  expect_equal(fit.lam4$beta[,1], lam4$beta[,1], tol=1e-4)
  expect_equal(edges4$edges$coef, expected=edgescoeff, tol=1e-3)
  })
