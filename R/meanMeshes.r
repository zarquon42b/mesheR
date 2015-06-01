#' average the vertices from a set of corresponding meshes
#'
#' average the vertices from a set of corresponding meshes
#' @param ... meshes of class mesh3d
#' @return a mesh with averaged vertices
#' @export
meanMeshes <- function(...) UseMethod("meanMeshes")

#' @export
meanMeshes.default <- function(...) {
    args <- list(...)
    ref <- meanMeshes(args)
    return(ref)
}

#' @export
meanMeshes.list <- function(...) {
    x <-  (...)
    mmmvert <- meshlist2array(x)
    mmmvertmean <- Morpho::arrMean3(mmmvert)
    ref <- x[[1]]
    ref$vb[1:3,] <- t(mmmvertmean)
    ref <- Rvcg::vcgUpdateNormals(ref)
    return(ref)
}


meshlist2array <- function (meshlist) {
    n <- length(meshlist)
    k <- ncol(meshlist[[1]]$vb)
    vertarr <- array(NA, dim = c(k, 3, n))
    for (i in 1:n) vertarr[, , i] <- vert2points(meshlist[[i]])
    dimnames(vertarr)[[3]] <- names(meshlist)
    if (is.null(names(meshlist))) 
        dimnames(vertarr)[[3]] <- paste("specimen", 1:n, sep = "_")
    return(vertarr)
}
