\name{mvregmed.edges}
\alias{mvregmed.edges}
\alias{plot.mvregmed.edges}
\title{
For an object of class mvregmed, create edges for a graph object that can be used
for plots, or for creating models input to lavaan function sem
}
\description{
Using the names of the alpha, beta, and delta matrices in the fitted
object, create directed edges 
}
\usage{
mvregmed.edges(fit, eps = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
    The fit as an object of class mvregmed. This can be output from
    either mvregmed.fit or mvregmed.grid.bestfit.
}
  \item{eps}{
    Threshold to determine whether any of alpha, beta, or delta
    parameters are close to zero to be rounded to zero.
}
}
\value{
  an object of class mvregmed.edges, which is a list with all.edge which
  is a data.frame with directed edges that are the names of x, mediator,
  y that are in the fitted model. Also returns data frames for alpha,
  beta, delta (each data frame containing vertex
  labels, row/col indices from which parameters were selected, and
  coefficient values), as well as names of x, mediator, and y.
}
\references{
 Schaid DJ, Dikilitas O, Sinnwell JP,  Kullo I. Penalized Mediation
 Models for Multivariate Data. Submitted, 2021.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
  \code{\link{mvregmed.fit}}
  \code{\link{plot.mvregmed.edges}}
}
\keyword{ model }
