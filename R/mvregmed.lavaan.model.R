mvregmed.lavaan.model <- function(fit.edge, fit.mvregmed)
  {
    mod <- NULL
    if(nrow(fit.edge$all.edge) < 1){
        stop("no edges in fit\n")
    } else {
        g.edge <- fit.edge$all.edge
    
        varx <- fit.mvregmed$varx
        varm <- fit.mvregmed$varm
        trait.name <- rownames(fit.mvregmed$vary)
    
        x.name <- rownames(varx)
        med.name <- rownames(varm)
    
        ## regression models
        dep.var <- unique(g.edge[, 2])
        for(i in 1:length(dep.var)){
            zed <- g.edge[,2] == dep.var[i]
            lhs <- dep.var[i]
            rhs <- g.edge[zed,1]
            rhs <- paste(rhs, collapse=" + ")
            new.mod <- paste(lhs, " ~ ", rhs, "\n", collapse="")
            mod <- paste(mod, new.mod, collapse="")
        }


        ## fixed med covar terms
        n.med <- length(med.name)
        for(i in 1:n.med){
            mod <- paste(mod, med.name[i], "  ~~ ",round(varm[i,i],3), " * ", med.name[i], " \n ", sep="")
        }
        for(i in 1:(n.med-1)){
            for(j in (i+1):n.med){
                mod <- paste(mod, med.name[i],"  ~~ ",round(varm[i,j],3), " * ", med.name[j], " \n ", sep="")
            }
        }

        ## fixed x covar terms
        n.x <- length(x.name)
        for(i in 1:n.x){
            mod <- paste(mod, x.name[i], "  ~~ ",round(varx[i,i],3), " * ", x.name[i], " \n ", sep="")
        }
        for(i in 1:(n.x-1)){
            for(j in (i+1):n.x){
                mod <- paste(mod, x.name[i],"  ~~ ",round(varx[i,j],3), " * ", x.name[j], " \n ", sep="")
            }
        }

        ## if there are any terms on lhs that are not on rhs and these terms
        ## are not traits (i.e., y), then the sem assumes that resid covar of these
        ## terms with traits  = 0.
        
        ulhs <- unique(fit.edge$all.edge[,2])
        urhs <- unique(fit.edge$all.edge[,1])
        which <- !(ulhs %in% urhs)
        ulhs <- ulhs[which]
        zed <- !(ulhs %in% trait.name)
        if(any(zed)){
            nm <- ulhs[zed]
            for(i in 1:length(nm)){
                for(j in 1:length(trait.name)){
                    mod <- paste(mod, nm[i],"  ~~ 0.0 * ", trait.name[j], " \n ", sep="")
                }
            }
        }
    }
    
    return(mod)
  }
