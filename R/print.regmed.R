print.regmed <- function(x, ...) {
 
  df <- data.frame(alpha=x$alpha, beta=x$beta)
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")  
  print(df, ...)
  invisible()
}
summary.regmed.grid <- function(object, ...) {
  print(object$grid.data, ...)
  invisible()
}
print.regmed.grid <- function(x, ...){
 
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$grid.data, ...)
  invisible()
}

