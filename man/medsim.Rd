\name{medsim}
\alias{medsim}
\alias{med}
\alias{x}
\alias{y}
\docType{data}
\title{
Simulated dataset for regmed package
}
\description{
Example data simulated from 2 response (y) variables, 10 exposure (x)
variables, and 200 mediators (med). The all variables are generated from
multivariate standard normal, with varying degrees of association
between exposure, mediators, and responses. 
}
\usage{data("medsim")}
\format{
  Three data frames with 100 observsations each:
  \describe{
    \item{\code{x}}{matrix with 10 columns of continuous exposure variables}
    \item{\code{y}}{matrix with 2 columns of continuous response variables}
    \item{\code{med}}{numeric matrix of 200 columns of simulated mediators between
      exposures and responses}
  }
}
\keyword{datasets}
