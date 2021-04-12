summary.lavaan <- function(object, ...){
    parm <- lavaan::parameterEstimates(object)
    zed <- !is.na(parm$z)
    parm <- parm[zed,]
    rownames(parm) <- NULL
    return(parm)
}
