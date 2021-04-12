\name{mvregmed.fit}
\alias{mvregmed.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate regularized mediation model 
}
\description{
  Fit regularized mediation model for a specified lambda penalty value.
  Structural equation models for analysis of multiple exposures (x),
  multiple mediators, and multiple outcome variables (y) are fit with a
  lasso (L1) penalaty on the model parameters.
  The model is x-[alpha] -> mediator-[beta] -> outcome, where
  alpha and beta are the parameters for the indirect effect of x on y,
  through the mediator. The model also allows a direct effect of x on y: x-[delta]->y.
}
\usage{
mvregmed.fit(x, mediator, y, lambda, x.std = TRUE, med.std = TRUE,
y.std = TRUE, max.outer = 5000, max.inner = 2, step.multiplier = 0.5,
print.iter = FALSE)
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
  
  \item{lambda}{
    lambda penalty parameter
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
  \item{max.outer}{
    maximum number of outer loop iterations. The outer loop cycles over several inner loops.
}
\item{max.inner}{
  maximum number of iterations for each inner loop. There is an inner loop
  for each paramemeter in the matrices alpha, beta, delta, and vary.
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
  An object of class mvregmed
}
\references{
  Schaid DJ, Dikilitas O, Sinnwell JP,  Kullo I. Penalized Mediation
 Models for Multivariate Data. Submitted, 2021.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
\code{\link{mvregmed.grid}}
}
\examples{
  data(medsim)
  mvfit <- mvregmed.fit(x, med[,1:10], y, lambda=.1)
  summary(mvfit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
