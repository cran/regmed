## Tests for mvregmed


context("Testing the mv-regmed")


data(medsim)

load("mvregmedfits.RData")

mvfit.grid <- mvregmed.grid(x, med[,1:20], y, lambda.vec=seq(.3, .04, by=-.01))


test_that("mvregmed_grid", {
  expect_equal(mvfit.grid$grid.data$df, grid20$grid.data$df)
  expect_equal(mvfit.grid$grid.data$bic, expected=grid20$grid.data$bic, tol=1e-2)
  })

mvfit.lam1 <- mvregmed.fit(x, med[,1:20], y, lambda=.1)

edges1 <- mvregmed.edges(mvfit.lam1)
deltas <- c(0.34227, .043031, -0.02459) # edges4$edges$coeff

#grid20 <- mvfit.grid
#lam1 <- mvfit.lam1
#save(grid20, lam1, file="mvregmedfits.RData")     


test_that("mvregmed_single", {
  expect_equal(mvfit.lam1$alpha[,1], lam1$alpha[,1], tol=1e-4)
  expect_equal(mvfit.lam1$beta[1,], lam1$beta[1,], tol=1e-4)
  expect_equal(edges1$delta.df$delta, expected=deltas, tol=1e-3)
  })
