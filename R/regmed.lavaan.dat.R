regmed.lavaan.dat <- function(x, mediator, y){
    dat.check <- regmed.dat.check(x, y, mediator)
    dat <- data.frame(cbind(scale(dat.check$x),scale(dat.check$mediator),scale(dat.check$y)))
    names(dat)[1] <- deparse(substitute(x))
    names(dat)[ncol(dat)] <- deparse(substitute(y))
    return(dat)
}
