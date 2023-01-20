\name{mvregmed.grid.data}
\alias{mvregmed.grid.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to summarize fits of models across a grid of lambda values
}
\description{
elper function to summarize fits of models across a grid of lambda values
}
\usage{
mvregmed.grid.data(fit.lst, lambda.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.lst}{
    A list of model fits over a grid of lambda values; length of list is
    length of vector of lambdas.
}
  \item{lambda.vec}{
    A vector of penalty lambda values/
}
}
\details{
  Create a data.frame of summmary information for each model fit in a grid.
}
\value{
  data.frame of summary information
}
\references{
  Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}

\seealso{
\code{\link{mvregmed.grid}}
}
\keyword{ model }
