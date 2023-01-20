mvregmed.lavaan.dat <- function(x, mediator, y, max.cor=0.99){
    dat.check <- mvregmed.dat.check(x, y, mediator, max.cor=max.cor)
    dat <- data.frame(cbind(scale(dat.check$x),scale(dat.check$mediator),scale(dat.check$y)))
    return(dat)
}
