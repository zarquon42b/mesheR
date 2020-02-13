#' Assign a color to a mesh globally
#'
#' Assign a color to a mesh globally
#'
#' @param mesh triangular mesh of class \code{mesh3d}
#' @param col color
#' @return
#' mesh with material adequately set.
#'
#' @examples
#' require(Rvcg)
#' data(humface)
#' humcol <- colorMesh(humface, "red")
#' @export colorMesh
colorMesh <- function(mesh, col) {
    if (!inherits(mesh,"mesh3d"))
        stop("please provide object of class mesh3d")

    col <- rgb(t(col2rgb(col[1])),maxColorValue = 255)
    material <- list()
    material$color <- rep(col,ncol(mesh$vb))
    mesh$material <- material
    invisible(mesh)
}

#' get the per-vertex colors of an object of class mesh3d
#'
#' get the per-vertex colors of an object of class mesh3d
#' @param x mesh3d
#' @return vector of hexadecimal rgb values
#' @export
getVertColor <- function(x) {
    if(is.null(x$material$color))
        stop("mesh has no per-vertex color")
    col <- 1:ncol(x$vb)
    tmp1 <- data.frame(it=as.vector(x$it))
    tmp1$rgb <- as.vector(x$material$color)
    tmp1 <- unique(tmp1)
    col[tmp1$it] <- tmp1$rgb
    return(col)
}

#' set the per-vertex colors of an object of class mesh3d
#'
#' set the per-vertex colors of an object of class mesh3d
#' @param x mesh3d
#' @param colorvec vector of lenght \code{ncol(x$vb)} containing HEX rgb values for vertex colors.
#' @return mesh with updated vertex colors
#' @export
setVertColor <- function(x,colorvec) {
    if (ncol(x$vb) != length(colorvec))
        stop("provide a color value for each vertex")
    colfun <- function(x){x <- colorvec[x];return(x)}
    x$material$color <- colorvec
    x$material$color[is.na(x$material$color)] <- #000000
    return(x)
}
