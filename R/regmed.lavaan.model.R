regmed.lavaan.model <- function(fit.edge, fit.regmed) {
    ## use edges object to construct structural equations for lavaan
    ## with lhs ~ rhs for all x-> med, and med-> y, and x->y
    ## that met the edges filtering applied to regmed.edges()
    mod <- NULL
    if(nrow(fit.edge$edges) < 1){
        stop("no edges in fit\n")
    }

    g.edge <- fit.edge$edges
    
    ## regression models
    dep.var <- unique(g.edge[, 2])
    med.keep <- NULL
    med.all <- rownames(fit.regmed$MedCov)
    for(i in 1:length(dep.var)){
        zed <- g.edge[,2] == dep.var[i]
        lhs <- dep.var[i]
        rhs <- g.edge[zed,1]
        med.keep <- unique(c(med.keep, med.all[med.all %in% c(lhs, rhs)]))          
        rhs <- paste(rhs, collapse=" + ")
        new.mod <- paste(lhs, " ~ ", rhs, "\n", collapse="")
        mod <- paste(mod, new.mod, collapse="")
    }
    
    
    ## fixed med covar terms
    
    ## determine which mediators selected in edges
    fit.edge$med.name %in% unique(as.vector(fit.edge$edges[, 1:2]))
    
    med.used <- fit.edge$med.name[fit.edge$med.name %in% unique(unlist(fit.edge$edges[, 1:2]))]
    
    select <- rownames(fit.regmed$MedCov) %in% med.used
    medcov <- fit.regmed$MedCov[select, select, drop=FALSE]
    n.used <- nrow(medcov)
    
    ## mediator with its variance
    for(i in 1:n.used){
        mod <- paste(mod,med.used[i], "  ~~ ",round(medcov[i,i],3), " * ", med.used[i], " \n ", sep="")
    }
    if(n.used > 1) {
        ## covariance of pairs of mediators
        for(i in 1:(n.used-1)){
            for(j in (i+1):n.used){
                mod <- paste(mod, med.used[i],"  ~~ ",round(medcov[i,j],3), " * ", med.used[j], " \n ", sep="")
            }
        }
    }

    ## if there are any terms on lhs that are not on rhs and these terms
    ## are not traits (i.e., y), then the sem assumes that resid covar of these
    ## terms with traits  = 0.
    
    ulhs <- unique(fit.edge$edges[,2])
    urhs <- unique(fit.edge$edges[,1])
    which <- !(ulhs %in% urhs)
    ulhs <- ulhs[which]
    zed <- !(ulhs %in% fit.edge$y.name)
    if(any(zed)){
        nm <- ulhs[zed]
        for(i in 1:length(nm)){
            mod <- paste(mod, nm[i],"  ~~ 0.0 * ", fit.edge$y.name, " \n ", sep="")
        }
    }

    return(mod)    
}

      
