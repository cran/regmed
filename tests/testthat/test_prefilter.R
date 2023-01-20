## Tests for prefilter functions for regmed and mvregmed

devel <- FALSE

context("Testing regmed.prefiler")

data(medsim)

if(devel) {
load("../../data/medsim.RData")
source("../../R/regmed.prefilter.R")  ## prefilters mediators
} 


set.seed(1000)
pfdat <- regmed.prefilter(x=x[,1], y=y[,1], mediator=med)
#sapply(pfdat, dim)
#colnames(pfdat1$mediator)
set.seed(111)
pfdatk5 <- regmed.prefilter(x=x[,1], y=y[,1], mediator=med, k=5)
#sapply(pfdatk5, dim)
#colnames(pfdatk5$mediator)
k5names <- paste0("med.", c(1,2,74,88,99))
test_that("regmed_prefilter", {
  expect_equal(ncol(pfdat$mediator), 50)
  expect_equal(colnames(pfdatk5$mediator), k5names)
  })

