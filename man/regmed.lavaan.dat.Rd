\name{regmed.lavaan.dat}
\alias{regmed.lavaan.dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set up data to input to lavaan sem
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Set up data to input to lavaan structural equation model (sem)
}
\usage{
regmed.lavaan.dat(x, mediator, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 vector of exposure variables
}
  \item{mediator}{
 matrix of mediator variables
}
  \item{y}{
 vector of outcome variable
}
}
\details{
 Use the function regmed.dat.check to standardize all variables and subset to subjects without missing data
}
\value{
A dataframe with updated x, mediator, and y
}
\references{
Schaid, DJ, Sinnwell JP. (2020) Penalized Models for Analysis of
Multiple Mediators. Genet Epidemiol 44:408-424.
}
\author{
 Dan Schaid, Jason Sinnwell
}
\keyword{ model }