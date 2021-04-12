\name{mvregmed.grid.bestfit}
\alias{mvregmed.grid.bestfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choose best fit model from a grid search based on minimum Bayesian
Information Criterion
}
\description{
Search over all models fit in a grid and choose model with min BIC as
best model.
}
\usage{
mvregmed.grid.bestfit(fit.grid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.grid}{
An object of class mvregmed.grid, output from function mvregmed.grid
}
}
\value{
An object of class mvregmed, a single best fitting model.
}
\references{
  Schaid DJ,Dikilitas O, Sinnwell JP,  Kullo I. Penalized Mediation
 Models for Multivariate Data. Submitted, 2021.
}
\author{
Daniel Schaid and Jason Sinnwell
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{mvregmed.grid}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
