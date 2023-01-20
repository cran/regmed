## Tests for dat.check functions for regmed and mvregmed

devel <- FALSE


context("Testing dat.check for regmed and mvregmed")

data(medsim)

if(devel) {
  load("../../data/medsim.RData")
  source("../../R/mvregmed.dat.check.R") ## scales data, gets rid of obs w/NA,
  source("../../R/regmed.dat.check.R")  ## scales and gets rid of obs w/NA
}


meddat1 <- regmed.dat.check(x=x[,1], y=y[,1], mediator=med, max.cor = 0.5)
sapply(meddat1, dim)


meddat2 <- mvregmed.dat.check(x=x, y=y, mediator=med, max.cor = 0.52)
sapply(meddat2, dim)


test_that("regmed_datcheck", {
  expect_equal(ncol(meddat1$mediator), 192)
  expect_equal(ncol(meddat2$mediator), 197)
  expect_equal(ncol(meddat2$x), 4)
  })
