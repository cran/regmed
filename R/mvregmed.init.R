mvregmed.init <- function(dat.obj,x.std=TRUE, med.std=TRUE, y.std=TRUE){
	
  ## Note that y should always be centered, so always do this, even if 
  ## input y was already standardized
  
  y <- scale(dat.obj$y,center=TRUE,scale=y.std)
  x <- scale(dat.obj$x,center=TRUE,scale=x.std)
  mediator <- scale(dat.obj$mediator,center=TRUE,scale=med.std) 

  nx <- ncol(dat.obj$x)
  nm <- ncol(dat.obj$mediator)
  ny <- ncol(dat.obj$y)

  dat <- cbind(x, mediator, y)
  sampcov <- var(dat)
 
  ## regress each med on all x to get residuals, using seemingly unrelated regression

  res <-apply(mediator,2,function(a,b) residuals(lm(a~0+b)),b=x)

  ## estimate penalized var matrix of residuals

  smat <- var(res)

  save.glasso <- glasso(smat, rho=.02, penalize.diagonal = FALSE)

  varm  <- save.glasso$w
  dimnames(varm)<-list(colnames(mediator),colnames(mediator))
    
  varx <- var(x)
  dimnames(varx) <- list(colnames(x), colnames(x))
    
  vary <- var(y)
  dimnames(vary) <- list(colnames(y), colnames(y))
    
  vary.step.size <- 0.05 * abs(min(vary))
  if(vary.step.size < 0.01) vary.step.size <- .01
    
  sampleSize <- nrow(dat)

  alpha <- matrix(0, nrow=nm, ncol=nx)
  dimnames(alpha) <- list(colnames(mediator), colnames(x))
    
  beta <- matrix(0, nrow= ny, ncol=nm)
  dimnames(beta) <- list(colnames(y), colnames(mediator))
   
  delta <- matrix(0, nrow=ny, ncol=nx)
  dimnames(delta) <- list(colnames(y), colnames(x))
 
  return(list(alpha=alpha,beta=beta,delta=delta,varx=varx, varm=varm,
              vary=vary,sampcov=sampcov,vary.step.size=vary.step.size,
              sampleSize=sampleSize))

}
