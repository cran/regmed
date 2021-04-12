summary.mvregmed <- function(object, eps = 1e-3, ...){

    df.tot <- sum(abs(object$alpha) > eps) + sum(abs(object$beta) > eps) + sum(abs(object$delta) > eps)
    if(df.tot == 0){
     warning("all alpha, beta, delta are zero")
     return (NULL)
    }

    
  object.edges <- mvregmed.edges(object, eps=eps)
  
  alpha.df <- object.edges$alpha.df
  alpha.df$row.index <- NULL
  alpha.df$col.index <- NULL
  
  beta.df <- object.edges$beta.df
  beta.df$row.index <- NULL
  beta.df$col.index <- NULL
  
  delta.df <- object.edges$delta.df
  delta.df$row.index <- NULL
  delta.df$col.index <- NULL
  cat("=== alpha parameter estimates ===\n")
  print(alpha.df, ...)
  cat("=== beta parameter estimates ===\n")
  print(beta.df, ...)
  cat("=== delta parameter estimates ===\n")
  print(delta.df, ...)

  invisible()
}
