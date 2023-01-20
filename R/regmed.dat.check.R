
regmed.dat.check <- function(x, y, mediator, max.cor=0.99, scale=TRUE){
  ## check input data for regmed (single x and y, many mediators)
    
  ###  mediator names ###
  mediator.names<-colnames(mediator)
  if(is.null(mediator.names)) mediator.names<-paste0("med.",1:ncol(mediator))
    
  ### check agreement of object dimensions ###
  if(any(c(length(x),length(y))!=nrow(mediator))) stop("dimensions of x,y and mediator do not agree")

  ## assure 0.01 < max.cor <=1
  if(max.cor < 0.01 | max.cor > 1.0) {
    warning(paste0("max.cor is outside allowed range: ", max.cor,
                   "\n  setting to defuault (0.99)\n"))
    max.cor <- 0.99
  }
    
  ### deal with missing data ###

  if(any(c(is.na(x),is.na(mediator),is.na(y)))) {

     if(options()$na.action=="na.pass") {

	stop("na.pass not allowed with missing data")

     } else {

	if(options()$na.action=="na.fail") {

		stop("missing data in x,y or mediator")

	} else {

		keep.obs <- !(is.na(y) | is.na(x) | rowSums(is.na(mediator)))
		x <- x[keep.obs]
		mediator <- mediator[keep.obs,]
		y <- y[keep.obs]

        }
      }
  }
    
  ## remove highly correlated mediators
  if(ncol(mediator) > 1 & max.cor < 1) {
      medmat <- abs(cor(mediator))
      idx.rc <- combinations(ncol(mediator),2)
      medvec <- medmat[idx.rc]
      idx.rc <- idx.rc[medvec > max.cor,,drop=FALSE]
      ## remove the second of the pair(s) that are correlated
      idx.rm <- sort(unique(idx.rc[,2]))
      if(length(idx.rm)>0) {
          warning("Note: removed due to high correlation, mediators: ",
              paste(colnames(mediator)[idx.rm], collapse=", "), "\n")
          mediator <- mediator[, !(1:ncol(mediator) %in% idx.rm), drop=FALSE]
          mediator.names <- mediator.names[!(1:ncol(mediator) %in% idx.rm)]
      }
  }
    
  return(list(x=x, y=y, mediator=mediator, mediator.names=mediator.names))

}
