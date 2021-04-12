\name{mvregmed.lavaan.dat}
\alias{mvregmed.lavaan.dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set up data to input to lavaan sem
}
\description{
Set up data to input to lavaan structural equation model (sem)
}
\usage{
mvregmed.lavaan.dat(x, mediator, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 matrix of exposure variables
}
  \item{mediator}{
 matrix of mediators
}
  \item{y}{
 matrix of outcome variables
}
}
\details{
Use the function regmed.dat.check to standardize all variables and subset to subjects without missing data
}
\value{
A dataframe with updated x, mediator, and y
}
\references{
  Schaid DJ,Dikilitas O, Sinnwell JP,  Kullo I. Penalized Mediation
 Models for Multivariate Data. Submitted, 2021.
}
\author{
Daniel Schaid and Jason Sinnwell
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
