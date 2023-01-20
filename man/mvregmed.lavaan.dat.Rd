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
mvregmed.lavaan.dat(x, mediator, y, max.cor=0.99)
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
\item{max.cor}{
  maximum correlation within mediators, so that fitting is more robust
}
}
\details{
Use the function regmed.dat.check to standardize all variables and subset to subjects without missing data
}
\value{
A dataframe with updated x, mediator, and y
}
\references{
  Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
