plot.mvregmed.grid <- function(x, xlab="lambda", ylab="BIC", pch="*", ...) {
  plot(x$grid.data$lambda, x$grid.data$bic,
       xlab=xlab, ylab=ylab, pch=pch, ...)
  invisible()
  
}
