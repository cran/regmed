\name{plot.mvregmed.grid}
\alias{plot.mvregmed.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot penalty parameter lambda versus BIC for model fits
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot penalty parameter lambda versus BIC for model fits
}
\usage{
\method{plot}{mvregmed.grid}(x, xlab="lambda", ylab="BIC", pch="*", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object created by mvregmed.grid
  }
  \item{xlab}{
    x axis label, by default set to 'lambda'
  }
  \item{ylab}{
    y axis label, by default set to 'BIC' for this plot
  }
  \item{pch}{
    plot character, by default set to a star (*)
  }
  \item{\dots}{
    optional plot arguments
  }
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
\keyword{model }
