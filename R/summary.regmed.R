summary.regmed <- function(object, ...){
  ## input:
  ## object is output from regmed.grid.best (i.e, a single fitted model) or 
  ## output from regmed.fit
   
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")  

  n.pts <- length(object$alpha)
   
  df <- data.frame(alpha=object$alpha, beta=object$beta,
                   alpha.beta=object$alpha * object$beta)
  colnames(df) <- c("alpha","beta","alpha*beta")
  print(df, ...)

  s.ab <- sum(object$alpha * object$beta)

  cat("\n")
  cat("sum of alpha*beta = ", s.ab, "\n")
  cat("delta = ", object$delta, "\n")
  cat("sum of delta + alpha*beta = ", object$delta + s.ab, "\n")
  cat("var(x) = ", object$var.x, "\n")
  cat("var(y) = ", object$var.y, "\n")
  invisible(df)
}
