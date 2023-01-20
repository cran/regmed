\name{summary.mvregmed}
\alias{summary.mvregmed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of mvregmed object
}
\description{
Summary of non-zero parameter estimates. Optional epsilon (eps) parameter controls
rounding to 0.
}
\usage{
\method{summary}{mvregmed}(object, eps=1e-3, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    mvregmed object returned from mvregmed.fit or mvregmed.grid.bestfit
  }
  \item{eps}{
    parameters smaller than epsilon (eps) are rounded to zero. 
  }
  \item{\dots}{
    optional arguments 
  }
}
\value{
Nothing is returned
}
\references{
  Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}

\seealso{
  \code{\link{mvregmed.fit}}
  \code{\link{mvregmed.grid.bestfit}}
}
\keyword{ model }
