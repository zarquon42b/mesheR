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
    material$color <- matrix(col, nrow(mesh$it), ncol(mesh$it))
    mesh$material <- material
    invisible(mesh)
}
