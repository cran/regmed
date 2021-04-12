## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE,
    fig.width=7, fig.height=6,
   tidy.opts=list(width.cutoff=100), tidy=TRUE, comment=NA)

## set lib path in R (below), 
#.libPaths(c("/projects/bsi/pharmacogenetics/s212953.SchaidMethods/Regmed/Version.2/Rlib362", .libPaths()))
## or with or with R_LIBS env variable
## linux: export R_LIBS=/projects/bsi/pharmacogenetics/s212953.SchaidMethods/Regmed/Version.2/Rlib362
library(igraph)
library(lavaan)
library(glasso)
library(regmed)

## ---- medsim------------------------------------------------------------------
data(medsim)


## ---- filtermediate-----------------------------------------------------------
dat.filter <- regmed.prefilter(x[,1], med, y[,1], k=10)

## ----set_grid-----------------------------------------------------------------

lambda.grid <- seq(from=.4, to=.01, by=-.01)

x1 <- dat.filter$x
y1 <- dat.filter$y
med <- dat.filter$mediator

fit.grid <- regmed.grid(x1, med, y1, lambda.grid, frac.lasso=.8)


## ---- methodsgrid-------------------------------------------------------------
plot.regmed.grid(fit.grid)

## ---- fitbest-----------------------------------------------------------------

fit.best <- regmed.grid.bestfit(fit.grid)
summary(fit.best)


## ---- singlefit---------------------------------------------------------------
fit.single <- regmed.fit(x1, med, y1, lambda=0.3, frac.lasso=.8)
summary(fit.single)

## ----edges--------------------------------------------------------------------

edges.med <- regmed.edges(fit.best, type="mediators")
plot.regmed.edges(edges.med)

edges.any <- regmed.edges(fit.best, type="any")
plot.regmed.edges(edges.any)


## ---- setuplavaan-------------------------------------------------------------

mod.best <- regmed.lavaan.model(edges.med, fit.best)

mod.any <- regmed.lavaan.model(edges.any, fit.best)

dat <- regmed.lavaan.dat(x1, med, y1)


## ---- fitlavaan---------------------------------------------------------------

fit.lav.med <- sem(model=mod.best, data=dat)
summary.lavaan(fit.lav.med)

fit.lav.any <- sem(model=mod.any, data=dat)
summary.lavaan(fit.lav.any)


## ----mvfit_grid---------------------------------------------------------------
fit.grid <- mvregmed.grid(x, med, y, lambda.grid)
summary(fit.grid)

## ----plot_mvgrid--------------------------------------------------------------
plot(fit.grid)

## ----mvfit_best---------------------------------------------------------------
mvfit.best <- mvregmed.grid.bestfit(fit.grid)


## ----mvfit_summary------------------------------------------------------------
summary(mvfit.best)

## ----fit_mvregmed-------------------------------------------------------------
mvfit <- mvregmed.fit(x, med, y, lambda=0.5)
summary(mvfit)

## ----plot_mvedges1------------------------------------------------------------

## get vertices and edges of graph 
mvfit.edges <- mvregmed.edges(mvfit.best, eps=1e-2)

plot.mvregmed.edges(mvfit.edges,
     x.color="palegreen",y.color="palevioletred",med.color="skyblue",
     v.size=20, seed=3)

## ----plot_mvedges2------------------------------------------------------------

gr <- regmed:::mvregmed.graph.attributes(mvfit.edges, x.color ="palegreen",
                                      y.color="palevioletred", 
                                      med.color="skyblue", v.size=30)
set.seed(3);
plot(gr$gr, vertex.size=gr$vsize, vertex.color=gr$vcol,
          vertex.label.font=1,vertex.label.color="black", 
          vertex.label.cex=.5,
          edge.arrow.mode=">",edge.arrow.size=.3)


## ----lavaanmv-----------------------------------------------------------------
mvmod <- mvregmed.lavaan.model(mvfit.edges, mvfit.best)
mvdat <- mvregmed.lavaan.dat(x, med, y)

mvfit.lavaan <- sem(model=mvmod, data=mvdat)

summary.lavaan(mvfit.lavaan)

