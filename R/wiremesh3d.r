#' create a shaded wiremesh
#'
#' create a shaded wiremesh
#' @param x mesh3d
#' @param col color of shaded mesh
#' @param specular parameter from rgl
#' @param ... additional parameters passed to shade3d from package rgl.
#' @return a vector containing the shape ids
#' @importFrom rgl wire3d shade3d
#' @rdname wiremesh3d
#' @export
wiremesh3d <- function(x,col="white",specular=1,...)UseMethod("wiremesh3d")

#' @rdname wiremesh3d
#' @export
wiremesh3d.mesh3d <- function(x,col="white",specular=1,...) {
    id1 <- wire3d(x)
    id2 <- shade3d(x,col=col,specular=specular,...)
    invisible( c(id1,id2))
}
