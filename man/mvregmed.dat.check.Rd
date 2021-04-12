\name{mvregmed.dat.check}
\alias{mvregmed.dat.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to check x, y, mediator for input to mvregmed functions 
}
\description{
Assure that x, y, mediator are matrices, check column names and fill in
if NULL, and reduce x, y, mediator so no missing values among all three matrices.
}
\usage{
mvregmed.dat.check(x, y, mediator)
}

\arguments{
  \item{x}{
    matrix with columns representing "exposure" variable (sometimes called instrumental variable)
}
  \item{y}{
    matrix with columns representing outcome variables
}
  \item{mediator}{
    matrix with columns representing mediator variables
}
}
\value{
list with updated x, y, mediator
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
\keyword{ models }

