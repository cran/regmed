\name{mvregmed.lavaan.model}
\alias{mvregmed.lavaan.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Setup a model for input to lavaan
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
   Set up a model statement (string formula) from mvregmed object for input to lavaan
}
\usage{
mvregmed.lavaan.model(fit.edge, fit.mvregmed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.edge}{
  Output from mvregmed.edges
}
  \item{fit.mvregmed}{
    A mvregmed object, either from mvregmed.fit or from mvregmed.grid.bestfit
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Loop through all relationships determined important from mrregmed edges
  object, and create model statement for lavaan, while also specifying
  covariances pre-estimated by mvregmed. See vignette for examples. The
  summary method for lavaan supersedes the summary from the lavaan
  package by simplifying the output to only return the coefficient
  table, as the covariance estimates were fixd from mvregmed.
}
\value{
Text string to define a model as input to sem
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
