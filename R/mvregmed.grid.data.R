mvregmed.grid.data <- function(fit.lst,lambda.vec){

    out<-data.frame(lambda=lambda.vec,
                    do.call(rbind,lapply(fit.lst,function(z) as.data.frame(z[c("converge","iter","df","df.alpha",
                                                                               "df.beta", "df.delta", "df.vary", "bic")]))))

    return(out)	

}
