mvregmed.dat.check <- function(x, y, mediator, max.cor=0.99){

    ## assure that x, y, mediator are matrices
    if(!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if(!is.matrix(y)) {
        y <- as.matrix(y)
    }
    if(!is.matrix(mediator)) {
        mediator <- as.matrix(mediator)
    }

    ## assure 0 < max.cor <=1
    if(max.cor < 0.01 | max.cor > 1.0) {
        warning(paste0("max.cor is outside allowed range: ", max.cor,
                       "\n",   "setting to default (0.99)\n"))
    }

    ##  mediator names ###
    mediator.names <-colnames(mediator)
    if(is.null(mediator.names)){
        mediator.names<-paste0("med.",1:ncol(mediator))
        colnames(mediator) <- mediator.names
    }
    
  ##  x names ###
  x.names <- colnames(x)
  if(is.null(x.names)){
      x.names<-paste0("x.",1:ncol(x))
      colnames(x) <- x.names
  }
    
  ##  y names ###
  y.names <- colnames(y)
  if(is.null(y.names)){
      y.names<-paste0("y.",1:ncol(y))
      colnames(y) <- y.names
  }
    
    
  ## check agreement of object dimensions ###
  if(any(c(nrow(x),nrow(y))!= nrow(mediator))){
      stop("row dim of x, y and mediator do not agree")
  }
    
  ## deal with missing data ###
    
  if(any(c(is.na(x),is.na(mediator),is.na(y)))) {
        
      if(options()$na.action=="na.pass") {
            
          stop("na.pass not allowed with missing data")
          
      } else {
          
          if(options()$na.action=="na.fail"){
              
              stop("missing data in x,y or mediator")
              
          } else {
              
              keep.obs <- !(rowSums(is.na(x))| rowSums(is.na(y))| rowSums(is.na(mediator)))
              
              x <- x[keep.obs,]
              mediator<-mediator[keep.obs,]
              y <- y[keep.obs,]
              
          }
      }
  }
    
  ## remove highly correlated (cor > max.cor) x exposures
  if(ncol(x) > 1) {
      cmat <- abs(cor(x))
      idx.rc <-  combinations(ncol(x),2)
      cvec <- cmat[idx.rc]
      idx.rc <- idx.rc[cvec > max.cor,,drop=FALSE]
      idx.rm <- sort(unique(idx.rc[,2]))
      if(length(idx.rm)>0) {
          warning("Note: removed due to high correlation, exposures (x): ",
                  paste(colnames(x)[idx.rm], collapse=", "), "\n")  
          x <- x[, !(1:ncol(x) %in% idx.rm), drop=FALSE]
      }
  }
    
  ## remove highly correlated mediators
  if(ncol(mediator) > 1 & max.cor < 1) {
      medmat <- abs(cor(mediator))
      idx.rc <- combinations(ncol(mediator),2)
      medvec <- medmat[idx.rc]
      idx.rc <- idx.rc[medvec > max.cor,,drop=FALSE]
      idx.rm <- sort(unique(idx.rc[,2]))
      if(length(idx.rm)>0) {
          warning("Note: removed due to high correlation, mediators: ",
              paste(colnames(mediator)[idx.rm], collapse=", "), "\n")
          mediator <- mediator[, !(1:ncol(mediator) %in% idx.rm), drop=FALSE]
      }
  }
    
  ## remove highly correlated y outcomes
  if(ncol(y) > 1) {
      cmat <- abs(cor(y))
      idx.rc <- combinations(ncol(y),2)
      cvec <- cmat[idx.rc]
      idx.rc <- idx.rc[cvec > max.cor,,drop=FALSE]
      idx.rm <- sort(unique(idx.rc[,2]))
      if(length(idx.rm)>0) {
          warning("Note: removed due to high correlation, outcomes (y): ",
              paste(colnames(y)[idx.rm], collapse=", "),
              "\n")  
          y <- y[, !(1:ncol(y) %in% idx.rm), drop=FALSE]
      }
  }
    
  return(list(x= x, y=y, mediator=mediator))

}
