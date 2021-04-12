
regmed.edges <- function(fit, type="mediators",  eps = 1e-3) {

    if(type != "mediators" && type != "any"){
        stop("unknown type of edges to select")
    }
    
    ## create directed edges for graph: vertex1 -> vertex2

    alpha <- fit$alpha
    beta <- fit$beta
    delta <- fit$delta

    df <- NULL
    
    if(type=="mediators"){
        df <- cbind(alpha, beta, alpha*beta)
        keep <- abs(df[,3]) > eps
        med.keep <- rownames(df[keep,, drop=FALSE])
        if(sum(keep) > 0) {
            df1 <- data.frame(Vertex1=colnames(alpha), Vertex2=med.keep, coef=alpha[keep])
            df2 <- data.frame(Vertex1=med.keep, Vertex2 = colnames(beta), coef=beta[keep])
            df <- rbind(df1, df2)
            if(abs(delta) > eps) {
                df3 <- data.frame(Vertex1 = colnames(alpha), Vertex2=colnames(beta),
                                  coef=delta)
                df <- rbind(df, df3)
                rownames(df) <- NULL
            }
        }
    }
        
    if(type == "any"){
             
        df1 <- data.frame(Vertex1=colnames(alpha), Vertex2=rownames(alpha), coef=alpha[,1])
        
        which <- abs(alpha) > eps
        if(any(which)){
            df <- df1[which,,drop=FALSE]
        }
        df2 <- data.frame(Vertex1=rownames(beta), Vertex2=colnames(beta), coef=beta[,1])
        which <- abs(beta) > eps
        if(any(which)){
            df <- rbind(df, df2 <- df2[which,,drop=FALSE])
            rownames(df)  <- NULL
        }
        if(abs(delta) > eps){
            df3 <- data.frame(Vertex1=colnames(alpha), Vertex2=colnames(beta), coef=delta)
            df <- rbind(df, df3)
            rownames(df) <- NULL
            }
        }

    y.name <- colnames(beta)
    med.name <- rownames(alpha)

        
    obj <- list(edges = df, y.name=y.name, med.name=med.name)
    class(obj) <- c("regmed.edges", "list")
    return(obj)
}

