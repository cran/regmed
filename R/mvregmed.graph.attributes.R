mvregmed.graph.attributes <- function(fit.edges, x.color ="palegreen",
                                      y.color="palevioletred", 
                                      med.color="skyblue", v.size=30){
  ## setup graph attributes
    edges <- as.matrix(fit.edges$all.edge)
    gr  <- graph_from_edgelist(edges, directed=TRUE)
    vname  <- vertex_attr(gr )$name
    vsize  <- rep(v.size, length(vname ))
    vcol  <- rep(x.color, length(vname))
    vcol[vname %in% fit.edges$y.name] <- y.color
    vcol[vname%in% fit.edges$med.name] <- med.color
    return(list(gr=gr, vname=vname, vsize=vsize, vcol=vcol))
 }
