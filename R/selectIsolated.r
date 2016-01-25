#' interactively select isolated pieces from a mesh
#'
#' interactively select isolated pieces from a mesh
#' @param mesh triangular mesh of class mesh3d
#' @param maxpiece integer: the n-largest (number of vertices) pieces to consider
#' @return returns a mesh with all selected parts merged and colored
#' @export
selectIsolated <- function(mesh,maxpiece=10) {
    pieces <- vcgIsolated(mesh,split=TRUE)
    ll <- length(pieces)
    nverts <- sapply(pieces,function(x) x <- ncol(x$vb))
    pieces <- pieces[order(nverts,decreasing = T)]
    pieces <- pieces[1:min(maxpiece,ll)]
    ll <- length(pieces)
    cols <- colorRampPalette(c("red","green","blue"))
    cols <- cols(ll)
    cat(paste0("there are ",ll," pieces\n"))
    open3d()
    out <- list()
    for (i in 1:ll) {
        wire3d(pieces[[i]],col=cols[i])
        answer <- readline("select (y/N)")
        if (answer %in% c("y","Y")) {
            out <- append(out,list(colorMesh(pieces[[i]],cols[i])))
        }
        else
            rgl.pop()
    }
    out <- mergeMeshes(out)
    return(out)
}
