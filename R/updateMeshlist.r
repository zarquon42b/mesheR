#' replace vertices in a meshlist with new values from an array
#'
#' replace vertices in a meshlist with new values from an array
#' @param x list containing triangular meshes of class mesh3d
#' @param array k x 3 x n array containing replacement vertices for x; the i-th matrix will be inserted in the i-th mesh
#' @return list with updated meshes
#' @export
updateMeshlist <- function(x,array) {
    if (length(x) != dim(array)[3])
        stop("array must contain vertices for each mesh")
    for (i in 1:length(x)) {
        x[[i]] <- updateVertices(x[[i]],array[,,i])
    }
    return(x)
}

#' replace vertices of a mesh
#'
#' replace vertices of a mesh
#' @param x mesh of class mesh3d
#' @param matrix with nrow(matrix) == amount of vertices in \code{x}
#' @return returns updated mesh
#' @export
updateVertices <- function(x,matrix) {
    if (ncol(x$vb) != nrow(matrix))
        stop("amount of new vertices must equal the amount of old vertices")
    x$vb[1:3,] <- t(matrix)
    x <- vcgUpdateNormals(x)
    return(x)
}
    
