\name{regmed.lavaan.model}
\alias{regmed.lavaan.model}
\alias{summary.lavaan}
\title{
Create a lavaan model
}
\description{
Create a lavaan model for input to lavan::sem()
}
\usage{
regmed.lavaan.model(fit.edge, fit.regmed)
}
\arguments{
  \item{fit.edge}{
    object created by regmed.edges function
  }
  \item{fit.regmed}{
    object created by regmed.fit from which the fit.edge object was made
  }
}
\value{
a character string that descibes the mediation model in format of lavaan model
}
\details{
  The fit.regmed object is needed for the fixed covariance estimates to
  be put into the model statement. The summary method supersedes the
  default summary from the lavaan package to only return the table of
  coefficients, as the covariances are fixed from regmed.fit.
}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\keyword{models }

