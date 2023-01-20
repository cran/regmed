## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE,
    fig.width=7, fig.height=6,
   tidy.opts=list(width.cutoff=100), tidy=TRUE, comment=NA)

require(regmed)


## ---- medsim------------------------------------------------------------------
data(medsim)

## ---- filtermediate-----------------------------------------------------------
dat.filter <- regmed.prefilter(x[,1], med, y[,1], k=10)
colnames(dat.filter$mediator)

## ----set_grid-----------------------------------------------------------------

lambda.grid <- seq(from=.4, to=.01, by=-.05)

x1 <- dat.filter$x
y1 <- dat.filter$y
med <- dat.filter$mediator

fit.grid <- regmed.grid(x1, med, y1, lambda.grid, frac.lasso=.8)


## ---- methodsgrid-------------------------------------------------------------
print(fit.grid)
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

dat <- regmed.lavaan.dat(x1, med, y1)

mod.best <- regmed.lavaan.model(edges.med, fit.best)

mod.any <- regmed.lavaan.model(edges.any, fit.best)


## ---- fitlavaan---------------------------------------------------------------

fit.lav.med <- lavaan:::sem(model=mod.best, data=dat)
summary.lavaan(fit.lav.med)

fit.lav.any <- lavaan:::sem(model=mod.any, data=dat)
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
mvedges <- mvregmed.edges(mvfit.best)
plot(mvedges)

## ----fit_mvregmed-------------------------------------------------------------
mvfit <- mvregmed.fit(x, med, y, lambda=0.4)
summary.mvregmed(mvfit, eps=0.01)

## ----plot_mvedges1------------------------------------------------------------

## get vertices and edges of graph 
mvfit.edges <- mvregmed.edges(mvfit.best, eps=5e-2)

plot.mvregmed.edges(mvfit.edges,
     x.color="palegreen",y.color="palevioletred",med.color="skyblue",
     v.size=20, seed=3)

## ----plot_mvedges2------------------------------------------------------------

gr <- regmed:::mvregmed.graph.attributes(mvfit.edges, x.color ="limegreen",
                                      y.color="palevioletred", 
                                      med.color="royalblue", v.size=40)
set.seed(3);
plot(gr$gr, vertex.size=gr$vsize, vertex.color=gr$vcol,
          vertex.label.font=2,vertex.label.color="black", 
          vertex.label.cex=.5,
          edge.arrow.mode=">",edge.arrow.size=.3)


## ----lavaanmv-----------------------------------------------------------------
mvmod <- mvregmed.lavaan.model(mvfit.edges, mvfit.best)
mvdat <- mvregmed.lavaan.dat(x, med, y)

mvfit.lavaan <- lavaan:::sem(model=mvmod, data=mvdat)

summary.lavaan(mvfit.lavaan)

