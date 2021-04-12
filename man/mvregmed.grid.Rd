\name{mvregmed.grid}
\alias{mvregmed.grid}
\alias{summary.mvregmed.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a grid of mvregmed models over a vector of lambda penalty parameters
}
\description{
For each lambda in an input vector of values, fit a penalized mvregmed
model
}
\usage{
mvregmed.grid(x, mediator, y, lambda.vec, max.outer = 5000,
max.inner = 2, x.std = TRUE, med.std = TRUE, y.std = TRUE,
step.multiplier = 0.5, print.iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    matrix with columns representing "exposure" variable (sometimes called instrumental variable)
  }
  \item{mediator}{
    matrix with columns representing mediator variables
    }
  \item{y}{
    matrix with columns representing outcome variables
  }
  \item{lambda.vec}{
   Vector of values of penalty parameter lambda's
}
 \item{max.outer}{
    maximum number of outer loop iterations. The outer loop cycles over several inner loops.
}
\item{max.inner}{
  maximum number of iterations for each inner loop. There is an inner loop
  for each paramemeter in the matrices alpha, beta, delta, and vary.
}
  \item{x.std}{
    logical (TRUE/FALSE) whether to standardize x by dividing by
    standard devation of x. Note that each column of x will be centered on its mean.
}
  \item{med.std}{
    logical (TRUE/FALSE) whether to standardize mediator by dividing by
    standard devation of mediator. Note that each column of mediator will be centered on its mean.
}
\item{y.std}{
  logical (TRUE/FALSE) whether to standardize y by dividing by
    standard devation of y. Note that each column of y will be centered on its mean.
}
\item{step.multiplier}{
   In inner loop, the step size is shrunk by the step.multiplier to
    assure that step size is not too large. Generally, the default of
    0.5 works well.
}
\item{print.iter}{
    print iteration number during fitting routine
  } 
}
\value{
An object of class mvregmed.grid
}
\references{
  Schaid DJ,Dikilitas O, Sinnwell JP,  Kullo I. Penalized Mediation
 Models for Multivariate Data. Submitted, 2021.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
\code{\link{mvregmed.fit}}
}
\examples{
  data(medsim)
  mvfit.grid <- mvregmed.grid(x, med[,1:10], y, lambda.vec=seq(.3, .04, by=-.01))
  summary(mvfit.grid)
  ## plot(mvfit.grid)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
