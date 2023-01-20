## Tests for regmed edges


context("Testing regmed edges")

load("regmedfits.RData")
load("mvregmedfits.RData")


rmedges.med <- regmed.edges(lam4)
#rmedges.med$edges
rmv1 <- c("filter5$x","med.1","filter5$x")

mvrmedges <- mvregmed.edges(lam1,  eps=0.001) 
#mvrmedges$all.edge
#mvrmedges$alpha.df
mvv1 <- c("x.1","x.1","x.3","x.1","x.8","x.10","med.1","med.10","med.2")

test_that("edges", {
    expect_equal(rmedges.med$edges[,"Vertex1"], expected=rmv1)
    expect_equal(mvrmedges$all.edge[,"vertex1"], expected=mvv1)
})

