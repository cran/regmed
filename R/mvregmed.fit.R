mvregmed.fit <-function(x, mediator, y, lambda, 
                        x.std =TRUE, med.std=TRUE, y.std=TRUE, max.outer=5000, max.inner=2,
                        step.multiplier = 0.5,
                        print.iter=FALSE){
  
  zed<-match.call()
 
  ## ----------input parameters
  ## x:  matrix of exposure variable
  ## mediator: matrix of possible mediators
  ## y: matrix of outcomes
  ## lambda: a single penalty parmeter
  ## max.outer: maximum number of iterations in outer loop of optimization steps
  ## max.inner: maximum number of iterations in inner loop of optimization steps
  ## x.std: if TRUE, standardize x before analyses (center and scale by standard deviation)
  ## med.std: if TRUE, standardize mediator before analyses
  ## print.iter: if TRUE, print when each iteration of optimization is conducted (verbose output)

  ### check bounded parameters  
  if(lambda < 0) stop("invalid lambda, must be >= 0")
  if(length(max.outer)!=1) stop("invalid max.outer, must be scalar")
  if(max.outer <= 0) stop("invalid max.outer, must be > 0")
  if(length(max.inner)!=1) stop("invalid max.inner, must be scalar")
  if(max.inner <= 0) stop("invalid max.inner, must be > 0")


  ### check x,y and mediator ###

  checked.dat <- mvregmed.dat.check(x=x,y=y,mediator=mediator)

  ### scale and center x,y and mediator as appropriate, initialize variables for rcpp_regmed ###

  inits <- mvregmed.init(dat.obj=checked.dat,x.std=x.std, med.std=med.std, y.std=y.std)
 
  ### fit regmed ###


  save <- rcpp_mvregmed(alpha=inits$alpha, beta=inits$beta, delta=inits$delta,
                       varx=inits$varx, varm=inits$varm, vary=inits$vary,
                       sampcov = inits$sampcov, sample_size= inits$sampleSize,
                       lambda = lambda, max_iter=max.outer, max_iter_inner=max.inner,
                       tol=1e-6, vary_step_size = inits$vary.step.size,
                       step_multiplier = step.multiplier,
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
         

 
  save$call <- zed
  
  class(save)<-"mvregmed"

  return(save)

}
