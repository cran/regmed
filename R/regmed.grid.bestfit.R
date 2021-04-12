regmed.grid.bestfit <- function(fit.grid){
  ## input:
  ## fit.grid is output of regmed.grid
    
  if(class(fit.grid) != "regmed.grid") {
      stop("input not regmed.grid class")
  }
    
  lambda <- fit.grid$grid.data$lambda
  bic <- fit.grid$grid.data$bic
  index.best <- (1:length(lambda))[bic == min(bic)]
  fit.best <- fit.grid$fit.list[index.best][[1]]  
  fit.best$MedCov <- fit.grid$MedCov
  class(fit.best) <- c("regmed", "list")
  return(fit.best) 
}
