\name{regmed.edges}
\alias{regmed.edges}
\alias{plot.regmed.edges}
\title{
For an object of class regmed, create edges for a graph object that can be used
for plots, or for creating models input to lavaan function sem
}
\description{
Using the names of the alpha, beta, and delta matrices in the fitted
object, create directed edges 
}
\usage{
regmed.edges(fit, type="mediators", eps = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
    The fit as an object of class mvregmed. This can be output from
    either mvregmed.fit or mvregmed.grid.bestfit.
  }
  \item{type}{
    Character string specifying whether to only keep edges for mediators
    that have a non-zero coefficient with exposure and response
    variables ("mediator") or all edges ("all").
  }
  \item{eps}{
    Threshold to determine whether any of alpha, beta, or delta
    parameters are close to zero to be rounded to zero.
}
}
\value{
  a list with class "regmed.edges" containing all.edge which is a
  data.framewith directed edges that are
  the names of x, mediator, y that are in the fitted model, with the
  coefficient for that edge. The plot method will plot the edges using
  igraph plotting options.
}
\references{
   Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
  \code{\link{regmed.fit}}
  \code{\link{plot.regmed.edges}}
}
\keyword{ model }
