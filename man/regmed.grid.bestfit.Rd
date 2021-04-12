\name{regmed.grid.bestfit}
\alias{regmed.grid.bestfit}

\title{
Find best fitting regmed model from regmed.grid object.
}
\description{
Find best fitting regmed model from regmed.grid object using minimum
BIC to select model.
}
\usage{
regmed.grid.bestfit(fit.grid)
}
\arguments{
  \item{fit.grid}{ a regmed.grid object
}
}
\value{
\item{ fit }{best fit regmed object based on minimum BIC}
\item{ grid }{row out of grid.data of regmed.grid object corresponding to best fit}
}
\references{
Schaid, DJ, Sinnwell JP. (2020) Penalized Models for Analysis of
Multiple Mediators. Genet Epidemiol 44:408-424.
}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\seealso{
\code{\link{regmed.grid}}
}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
