mvregmed.edges <- function(fit, eps = 1e-3) {

    ## create directed edges for graph: vertex1 -> vertex2

    all.edge <- NULL
    a.df <- b.df <- d.df <- NULL
    
    mat <- cbind(row(fit$alpha)[abs(fit$alpha) > eps],
                 col(fit$alpha)[abs(fit$alpha) > eps])

    df.tot <- sum(abs(fit$alpha) > eps) + sum(abs(fit$beta) > eps) + sum(abs(fit$delta) > eps)
    if(df.tot == 0){
     warning("no edges, all alpha, beta, delta are zero")
     return()
    }
    
    
    if(nrow(mat) > 0){
        ord <- order(mat[,1], mat[,2])
        mat <- mat[ord,]

        a.df <- data.frame(mediator=rownames(fit$alpha)[mat[,1]],
                           x=colnames(fit$alpha)[mat[,2]],
                           row.index=mat[,1],
                           col.index=mat[,2],
                           alpha=fit$alpha[mat])
        all.edge <- rbind(all.edge, as.matrix(a.df[, c("x","mediator")]))
        }
    
    mat <- cbind(row(fit$delta)[abs(fit$delta) > eps],
                 col(fit$delta)[abs(fit$delta) > eps])
    
    if(nrow(mat) > 0){
        ord <- order(mat[,1], mat[,2])
        mat <- mat[ord,]
    
        d.df <- data.frame(y=rownames(fit$delta)[mat[,1]],
                           x=colnames(fit$delta)[mat[,2]],
                           row.index=mat[,1],
                           col.index=mat[,2],
                           delta=fit$delta[mat])
        all.edge <- rbind(all.edge, as.matrix(d.df[,c("x","y")]))
    }
    
    mat <- cbind(row(fit$beta)[abs(fit$beta) > eps],
                 col(fit$beta)[abs(fit$beta) > eps])
    if(nrow(mat) > 0){
        ord <- order(mat[,1], mat[,2])
        mat <- mat[ord,]
    
        b.df <- data.frame(y=rownames(fit$beta)[mat[,1]],
                           mediator=colnames(fit$beta)[mat[,2]],
                           row.index=mat[,1],
                           col.index=mat[,2],
                           beta=fit$beta[mat])
        all.edge <- rbind(all.edge,  as.matrix(b.df[,c("mediator","y")]))
    }

    if(is.null(all.edge)){
        warning("no edges created")
        return()
    }

       
    if(ncol(all.edge) == 2){
        colnames(all.edge) <- c("vertex1", "vertex2")
    }
    
    obj <- list(all.edge=all.edge,
                alpha.df=a.df,
                beta.df=b.df,
                delta.df=d.df,
                x.name = colnames(fit$alpha),
                med.name = rownames(fit$alpha),
                y.name = rownames(fit$delta))
    class(obj) <- c("mvregmed.edges", "list")
    return(obj)
  }

