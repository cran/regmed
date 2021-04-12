\name{plot.regmed.grid}
\alias{plot.regmed.grid}
\title{
Plots for regmed.grid object.
}
\description{
Creates 2 plots: (1) BIC vs. lambda, and (2) Coefficients Alpha and Beta of mediator vs. lambda.
}
\usage{
\method{plot}{regmed.grid}(x, as.log=FALSE, \dots)
}
\arguments{
  \item{x}{ regmed.grid object, returned by the \code{regmed.grid()} function }
  \item{as.log}{Logical; if TRUE,  plot lambda on the log scale}
  \item{\dots}{optional arguments for plot method}
}
\value{nothing is returned}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\seealso{
\code{\link{regmed.grid}},
}
\examples{
data(medsim)
fit.grid <- regmed.grid(x[,1], med[,1:10], y[,1],
  lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
summary(fit.grid)
}
