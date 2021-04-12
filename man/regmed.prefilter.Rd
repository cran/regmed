\name{regmed.prefilter}
\alias{regmed.prefilter}
\alias{regmed.dat.check}
\title{
  Prefilter to reduce the number of mediators for subsequent analyses
}
\description{
Use sure independence screening (Fan & Lv, 2008)to reduce the number of
potential mediators when the number of potential mediators is
large. This is based on ranking marginal correlations and then selecting
the highest ranked values such that the number of parameters is less
than the sample size. Because mediation depends on the two correlations,
cor(x,med) and cor(med, y) we rank the absolute values of their
products, |cor(x, med) * cor(med, y)|, and choose the highest k ranked
values to determine which potential mediators to include in penalized
mediation models. If k is not specified, the default value of k is n/2,
where n is the sample size, because each mediator results in two
parameters alpha and beta.
}
\usage{
regmed.prefilter(x, mediator, y, k = NULL, x.std = TRUE,
med.std = TRUE, y.std=TRUE)
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
\item{k}{
  Number of potential mediators to select. Default is n/2, where n is sample size.
}
\item{x.std}{ 
  logical (TRUE/FALSE) whether to standardize x by dividing by standard
  devation of x. Note that x will be centered on its mean.
}
\item{med.std}{
  logical (TRUE/FALSE) whether to standardize mediators by dividing each
  mediator by its standard deviation. Note that mediators will be centered on their means.
}
\item{y.std}{ 
  logical (TRUE/FALSE) whether to standardize y by dividing by standard
  devation of y. Note that y will be centered on its mean.
}
}
\value{
list with x, mediator, and y, after subsetting to no missing values,
applying x.std and med.std, and subsetting mediators to k top choices.
}
\references{
Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature space. J. R. Statist. Soc.B, 70, 849-911. 
Schaid, DJ, Sinnwell JP. (2020) Penalized Models for Analysis of
Multiple Mediators. Genet Epidemiol 44:408-424.
}
\author{
  Dan Schaid, Jason Sinnwell
}
\seealso{
\code{\link{regmed}}
}
\examples{
data(medsim)
dim(med)
filtered <- regmed.prefilter(x[,1], med, y[,1], k=10)
dim(filtered$med)
}
\keyword{ models }
