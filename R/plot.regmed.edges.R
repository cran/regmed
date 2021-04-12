
plot.regmed.edges <- function(x, x.color="palegreen",
          y.color="palevioletred", med.color="skyblue",
          vertex.label.color="black", v.size=30, seed=100, ...) {


    ## create graph attributes: color of vertices and size
    gredges <- as.matrix(x$edges[,1:2])
    gr  <- graph_from_edgelist(gredges, directed=TRUE)
  #  browser()
    vname  <- vertex_attr(gr)$name
    vsize  <- rep(v.size, length(vname ))
    vcol  <- rep(x.color, length(vname))
    vcol[vname %in% x$y.name] <- y.color
    vcol[vname%in% x$med.name] <- med.color


    ## igraph generates graphs based on random seed, so set seed
    ## to assure plot can be recreated. See documentation for igraph
    ## for further graphical parameters.

    set.seed(seed)
    plot(gr, vertex.size=vsize, vertex.color=vcol,
          vertex.label.font=1, vertex.label.color=vertex.label.color, 
          vertex.label.cex=.5, edge.arrow.mode=">", edge.arrow.size=.3)

    invisible()
}
