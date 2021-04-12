mvregmed.grid.update <- function(fit.obj, inits){

    inits$alpha <- fit.obj$alpha
    inits$beta  <- fit.obj$beta
    inits$delta <- fit.obj$delta
    inits$varx  <- fit.obj$varx
    inits$varm  <- fit.obj$varm
    inits$vary  <- fit.obj$vary

    ## avoid very small vary along diag
    dv <- diag(inits$vary)
    dv <- ifelse(dv < .5, .5, dv)
    diag(inits$vary) <- dv
    
   return(inits)

}
