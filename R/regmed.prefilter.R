## Pre-filter mediators  on N subjects to K < N/2, or user threshold

regmed.prefilter <- function(x, mediator, y, k=NULL, x.std=TRUE,
                             med.std=TRUE, y.std=TRUE, max.cor=0.99) {

    ## x exposure (numeric)
    ## mediator (numeric)
    ## y outcome/response variable (numeric)

    ## return list of x, y, mediator that is subset to k <= n/2 mediators,
    ## all subset to complete non-missings
    
    ## handle missings and too high of correlations
    checked.dat <- regmed.dat.check(x = x, y = y, mediator = mediator, max.cor=max.cor)

    y.std <- scale(checked.dat$y, center=TRUE, scale=y.std)
    x.std <- scale(checked.dat$x, center=TRUE, scale=x.std)
    mediator.std <- scale(checked.dat$mediator,center=TRUE, scale=med.std) 

    
    ## how many mediators
    if(is.null(k)) {       
        if(ncol(mediator.std) <= nrow(mediator.std)/2) {
            cat("filtering not needed (k <= nsubj/2)\n")
            return(list(x=checked.dat$x,
                        mediator=checked.dat$mediator,
                        y=checked.dat$y))
        } else {
            k <- ceiling(nrow(mediator.std)/2)
        }
    }
    
    ## to narrow down the mediators based on correlation to x and y
    r.x <- cor(x.std, mediator.std)
    r.y <- cor(y.std, mediator.std)

    r.xy <- abs(r.x*r.y)
    ##  choose those with higher value for r.xy, which is what rank does
    rank.meds <- rank(r.xy)
    idx.keep <- which(rank.meds > (ncol(mediator) - k))
    
    ## return list of x, y, mediator filtered
    return(list(x=checked.dat$x,
                mediator=checked.dat$mediator[,idx.keep],
                y=checked.dat$y))
}
    
