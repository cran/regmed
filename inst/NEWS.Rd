\name{NEWS}
\title{NEWS file for the regmed package}
\section{Release Notes for version 2.1}{
  \itemize{
    \item add max.cor for maximum correlation of mediators in regmed,
    and all x,y,mediators in mvregmed. Penalization and fitting struggle with highly correlated variables.
    \item fix bug in mvregmed.edges when down to one row for either  alpha, beta, delta
    \item remove procedural notes in regmed.grid and mvregmed.grid 
  }
}
\section{Release Notes for version 2.0}{
  \itemize{
    \item add mvregmed.grid and mvregmed.fit
    \item add functions to create edges of DAG from fitted object\
    \item add functions to plot the edges as a DAG
    \item functions to make running and summarizing lavaan runs
    \item updated vignette with mvregmed examples
    \item replace trim.best functionality with pulling fit with best bic lambda
  }
}
\section{Release Notes for version 1.1}{
  \itemize{
    \item Fix BIC calculation
    \item default lambda=0 for trim.best()
  }
} 
\section{Release Notes for version 1.0}{
  \itemize{
    \item Original release to CRAN
    \item Regularized Mediation analysis with plot methods for
    regmed.fit and regmed.grid S3 class objects.
    \item Vignette that demonstrates usage with four mediators.
  }
}    

