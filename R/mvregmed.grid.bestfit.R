mvregmed.grid.bestfit <- function(fit.grid){

    if(!("mvregmed.grid" %in%  class(fit.grid))){
        stop("object not mvregmed.grid class")
    }
    
    bic <- fit.grid$grid.data$bic
    index <- 1:length(bic)
    lambda <- fit.grid$grid.data$lambda
    zed <- bic == min(bic)
    index.select <- index[zed]
    
    ## if tied for bic, choose model with largest lambda
    if(sum(zed) > 1){
        index.select <- index.select[lambda[zed] == max(lambda[zed])]
    }
    
   fit.best <- fit.grid$fit.list[[index.select]]
   class(fit.best)<-"mvregmed"
   return(fit.best)
 }
