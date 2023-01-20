\name{mvregmed.init}
\alias{mvregmed.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to setup data and parameters for input to mvregmed.fit
and mvregmed.grid
}
\description{
Helper function to setup data and parameters for input to mvregmed.fit
and mvregmed.grid
}
\usage{
mvregmed.init(dat.obj, x.std = TRUE, med.std = TRUE, y.std = TRUE)
}
\arguments{
  \item{dat.obj}{
    A list that is output from mvregmed.dat.check that contains x, mediator, and y.
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

}
\details{
Center and scale (if declared) x, mediator and y. Then regress each
mediator on all x to create residuals that are used to create the
residual variance matrix for mediators. This variance matrix is
penalized by glasso to obtain a matrix of full rank. Variance matrices
for x and y variables are also created. Initial values of paramemeter
matrices alpha, beta, and delta are created (all intital values = 0).
}
\value{
A list of items used as input to model fitting.
}
\references{
  Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized
  mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
  \code{\link{mvregmed.fit}}
  \code{\link{mvregmed.grid}}
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
