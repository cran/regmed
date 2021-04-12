
plot.mvregmed.edges <- function(x, x.color="palegreen",
          y.color="palevioletred", med.color="skyblue",
          vertex.label.color="black", v.size=30, seed=100, ...) {


    ## create graph attributes: color of vertices and size
    gr <- mvregmed.graph.attributes(x, x.color ="palegreen",
                                     y.color="palevioletred", 
                                     med.color="skyblue", v.size=v.size, ...)

    ## igraph generates graphs based on random seed, so set seed
    ## to assure plot can be recreated. See documentation for igraph
    ## for further graphical parameters.

    set.seed(seed)
    plot(gr$gr, vertex.size=gr$vsize, vertex.color=gr$vcol,
          vertex.label.font=1, vertex.label.color=vertex.label.color, 
          vertex.label.cex=.5, edge.arrow.mode=">", edge.arrow.size=.3)

    invisible()
}
