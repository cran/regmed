\name{mvregmed.graph.attributes}
\alias{mvregmed.graph.attributes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setup attributes of graph object for plotting
}
\description{
Setup attributes of graph object for plotting directed acyclic graph
with graph attributes
}
\usage{
mvregmed.graph.attributes(fit.edges, x.color = "palegreen",
y.color = "palevioletred", med.color = "skyblue", v.size = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.edges}{
   A data.frame of edges with 1st column a vertex directed to the
   vertex in the 2nd column. This is all.edge from the list
   that is output from mvregmed.edges
}
  \item{x.color}{
    Color of vertices for x variables
}
  \item{y.color}{
    Color of vertices for y variables
}
  \item{med.color}{
    Color of vertices for mediators
}
  \item{v.size}{
    Size of vertices for plotting
}
}
\details{
  User can use this as template for taking advantage of more igraph attributes.
}
\value{
List with items 1) output from graph_from_edgelist (see igraph); 2) vertex names;
3)vertex size; 4) vertex colors
}
\references{
  Schaid DS, Dikilitas O, Sinnwell JP, Kullo I (2022). Penalized mediation models for multivariate data. Genet Epidemiol 46:32-50.
}
\author{
Daniel Schaid and Jason Sinnwell
}
\seealso{
  \code{\link{mvregmed.edges}}
  \code{\link{graph_from_edgelist}}
}
\keyword{ model}
