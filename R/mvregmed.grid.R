## this version breaks out of lambda loop if df > n

mvregmed.grid  <- function(x, mediator, y, lambda.vec,  max.outer=5000, max.inner=2,
                          x.std=TRUE, med.std=TRUE, y.std=TRUE, step.multiplier = 0.5,
                          print.iter = FALSE){
  
  zed <- match.call()

  
  ## ----------input parameters
  ## x:  vector of exposure variable
  ## mediator: matrix of possible mediators
  ## y: vector of outcome
  ## lambda.vec: vector of penalty parmeters (preferable in decreasing order for
  ##             using parameters estimates from one penalty output as warm start
  ##             initial values for next penalty parameter)
  ## frac.lasso: fraction of penalty assigned to lasso penalty (in contrast to
  ##             group penalty). Simulations suggest frac.lasso = 0.8 offers
  ##             adequate control of false-positive results without much loss in power
  ## max.outer: maximum number of iterations in outer loop of optimization steps
  ## max.inner: maximum number of iterations in inner loop of optimization steps
  ## x.std: if TRUE, standardize x before analyses (center and scale by standard deviation)
  ## med.std: if TRUE, standardize mediator before analyses
  ## print.iter: if TRUE, print when each iteration of optimization is conducted (verbose output)

  tol <- 1e-4 


  ### check bounded parameters  
  if(any(lambda.vec < 0)) stop("invalid lambda, lambda must be >= 0")
  if(length(max.outer)!=1) stop("invalid max.outer, must be scalar")
  if(max.outer <= 0) stop("invalid max.outer, must be > 0")
  if(length(max.inner)!=1) stop("invalid max.inner, must be scalar")
  if(max.inner <= 0) stop("invalid max.inner, must be > 0")

  ### check x,y and mediator ###

  checked.dat <- mvregmed.dat.check(x=x,y=y,mediator=mediator)

  ### scale and center x,y and mediator as appropriate, initialize variables for rcpp_regmed ###

  inits <- mvregmed.init(dat.obj=checked.dat, x.std=x.std, med.std=med.std, y.std=y.std)
   
  ### loop through lambda values ###

  fit.lst <- vector("list", length=length(lambda.vec))

  lambda.used <- NULL

    
  for(i in 1:length(lambda.vec)){
  
    cat("fitting grid lambda [",i,"] = ", lambda.vec[i], "\n")

      save <- rcpp_mvregmed(alpha=inits$alpha, beta=inits$beta, delta=inits$delta,
                       varx=inits$varx, varm=inits$varm, vary=inits$vary,
                       sampcov = inits$sampcov, sample_size= inits$sampleSize,
                       lambda = lambda.vec[i], max_iter=max.outer, max_iter_inner=max.inner,
                       tol=1e-6, step_multiplier = step.multiplier,
                       verbose=print.iter)
      
      ## add names to parameters
                        
        rownames(save$alpha) <- colnames(checked.dat$mediator)
        colnames(save$alpha) <- colnames(checked.dat$x)
        
        rownames(save$beta) <- colnames(checked.dat$y)
        colnames(save$beta) <- colnames(checked.dat$mediator)
        
        rownames(save$delta) <- colnames(checked.dat$y)
        colnames(save$delta) <- colnames(checked.dat$x)
        
        rownames(save$varx)  <- colnames(checked.dat$x)
        colnames(save$varx)  <- colnames(checked.dat$x)
        
        rownames(save$vary)  <- colnames(checked.dat$y)
        colnames(save$vary)  <- colnames(checked.dat$y)
        
        rownames(save$varm)  <- colnames(checked.dat$mediator)
        colnames(save$varm)  <- colnames(checked.dat$mediator)
         
    ### zero out alphas/betas/deltas if they are sufficiently close to 0, based on tol parameter ###    

    save$alpha <- ifelse(abs(save$alpha) < tol, 0.0, save$alpha)
    save$beta  <- ifelse(abs(save$beta)  < tol, 0.0, save$beta)
    save$delta  <- ifelse(abs(save$delta)  < tol, 0.0, save$delta)

      fit.lst[[i]] <- save

      ## update inits for next loop value of lambda: alpha, delta, beta, varx, varm, var

      inits <- mvregmed.grid.update(fit.obj=save, inits=inits)

      ## check on invalid alpha/beta/delta if not on last lambda value, to avoid
      ## using invalid parameters as warm starts for fit of next lambda.
      ## This allows the last lambda value to possibly result in invalid parameters

        lambda.used <- c(lambda.used, lambda.vec[i])
      
      if(i < length(lambda.vec)){
          if(any(is.nan(save$alpha)) | any(is.na(save$alpha)) ) stop("invalid alpha estimates")
          if(any(is.nan(save$beta))  | any(is.na(save$beta))  ) stop("invalid beta estimates")
          if(any(is.nan(save$delta)) | any(is.na(save$delta)) ) stop("invalid delta estimates")
      }

      ## break out of loop if number of estimated parameters (df) > sample size 
      if(save$df > nrow(x)){
          break
          }
  }
    

      ## pass lambda.used in case shorter than lambda.vec (if break out of loop)

    fit.lst <- fit.lst[1:length(lambda.used)]
  
  gridData <- mvregmed.grid.data(fit.lst=fit.lst,lambda.vec=lambda.used)

  
  out<-list(fit.list=fit.lst, grid.data=gridData, sample.size=inits$sampleSize, MedCov=inits$MedCov, call=zed)
  class(out)<-"mvregmed.grid"

  return(out)

}

summary.mvregmed.grid <- function(object, ...){
  print(object$grid.data, ...)
}
