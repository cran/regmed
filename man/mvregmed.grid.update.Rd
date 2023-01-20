\name{mvregmed.grid.update}
\alias{mvregmed.grid.update}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to update parameters in a grid search
}
\description{
After a model is fit with a specific lambda, use the output of the
fitted parameters as initial values for the next lambda value, thus
using warm starts at each successive lambda value.
}
\usage{
mvregmed.grid.update(fit.obj, inits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{
A fitted model of class mvregmed.
}
  \item{inits}{
    Initial values from mvregmed.init that are subsequently updated with
    new values from fit.obj
  }
}
\value{
  A list with the same components as output from mvregmed.init
}
\references{
  Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
\code{\link{mvregmed.init}}
}
\keyword{ model }
