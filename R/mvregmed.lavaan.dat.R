mvregmed.lavaan.dat <- function(x, mediator, y){
    dat.check <- mvregmed.dat.check(x, y, mediator)
    dat <- data.frame(cbind(scale(dat.check$x),scale(dat.check$mediator),scale(dat.check$y)))
    return(dat)
}
