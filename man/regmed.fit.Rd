\name{regmed.fit}
\alias{regmed.fit}
\alias{print.regmed}
\alias{summary.regmed}
\title{
  Regularized Mediation model for a specified lambda penalty value.
}
\description{
Fit regularized mediation model for a specified lambda penalty value.
Structural equation models for analysis of multiple mediators
are extended by creating a sparse group lasso penalized model such that
the penalty considers the natural groupings of the pair of parameters
that determine mediation, as well as encourages sparseness of the model
parameters. The model is x-[alpha] -> mediator-[beta] -> outcome, where
alpha and beta are the parameters for the indirect effect of x on y,
through the mediator. The model also allows a direct effect of x on y: x-[delta] -> y.
}
\usage{
regmed.fit(x, mediator, y, lambda, frac.lasso=0.8, x.std=TRUE, med.std=TRUE,
max.outer=5000, max.inner=100, step.multiplier = 0.5, wt.delta = .5,
print.iter=FALSE, max.cor=0.99)
}
\arguments{
  \item{x}{
vector representing "exposure" variable (sometimes called instrumental variable)
}
  \item{mediator}{
matrix of mediators, rows are observations, columns are different mediators
}
  \item{y}{
vector representing outcome
}
  \item{lambda}{
lambda penalty parameter
}
  \item{frac.lasso}{
fraction of penalty (lambda) that is allocated to L1 penalty
(lasso). The remaining fraction, (1-frac.lasso) is allocated to
group-lasso penalty, where the group is the pair of parameters alpha and
beta that determine mediation (x [alpha] -> mediator -> [beta] y).
}
  \item{x.std}{ 
logical (TRUE/FALSE) whether to standardize x by dividing by standard
devation of x. Note that x will be centered on its mean.
}
  \item{med.std}{
logical (TRUE/FALSE) whether to standardize mediators by dividing each
mediator by its standard deviation. Note that mediators will be centered on their means.
}
  \item{max.outer}{
maximum number of outer loop iterations. The outer loop cycles over several inner loops.
} 
  \item{max.inner}{
maximum number of iterations for each inner loop. There is an inner loop
for each pair of alpha-beta parameters for each mediator, an inner loop 
for direct effect (delta), and inner loops for residual variances for x
and for y.
}
 \item{step.multiplier}{
a value between 0 and 1 for backtracking, to shrink step size. Value of
0.5 is typical default.
}
\item{wt.delta}{
  a weight >=0 for how much weight should be given to shrinking delta
  parameter, by penalaty lambda*wt.delta.
}

\item{print.iter}{
print iteration history during fitting routine
}
\item{max.cor}{
  maximum correlation within mediators, so that fitting is more robust
}
}
\value{
regmed object, with S3 methods available: plot, print, summary
}
\references{
Schaid, DJ, Sinnwell JP. (2020) Penalized Models for Analysis of
Multiple Mediators. Genet Epidemiol 44:408-424.
}
\author{
  Dan Schaid, Greg Jenkins, Jason Sinnwell
}

\seealso{
\code{\link{regmed.edges}}
\code{\link{summary.regmed}}
}
\examples{
  data(medsim)
  filter5 <- regmed.prefilter(x[,1], med, y[,1], k=5)
  fit.regmed <- regmed.fit(x[,1], med[,1:5], y[,1], lambda = 0.2, frac.lasso=.8)
  summary(fit.regmed)
}
\keyword{ models }
