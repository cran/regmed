mvregmed.dat.check <- function(x,y,mediator){

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

    ##  mediator names ###
    mediator.names <-colnames(mediator)
    if(is.null(mediator.names)){
        mediator.names<-paste0("med.",1:ncol(mediator))
        colnames(mediator) <- mediator.names
    }
    
  ##  x names ###
  x.names<-colnames(x)
  if(is.null(x.names)){
      x.names<-paste0("x.",1:ncol(x))
      colnames(x) <- x.names
      }
    
    ##  y names ###
    y.names<-colnames(y)
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
            
        }else{
            
            if(options()$na.action=="na.fail"){
                
		stop("missing data in x,y or mediator")
                
            }else{
                
		keep.obs <- !(rowSums(is.na(x))| rowSums(is.na(y))| rowSums(is.na(mediator)))
                
		x<-x[keep.obs,]
		mediator<-mediator[keep.obs,]
		y<-y[keep.obs,]
                
            }
        }
    }
    
    return(list(x=x,y=y,mediator=mediator))

}
