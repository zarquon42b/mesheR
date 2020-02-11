#' inflate a mesh along its normals
#' 
#' translate vertices of a triangular mesh along its normals
#' 
#' 
#' @param mesh triangular mesh of class "mesh3d"
#' @param offset distance to translate the vertices
#' @return returns modified mesh.
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' require(Morpho)
#' data(nose)
#' offset <- meshOffset(shortnose.mesh,3)
#' \dontrun{
#' require(rgl)
#' shade3d(shortnose.mesh,col=3)
#' wire3d(offset)
#' }
#' @export meshOffset 
meshOffset <- function(mesh,offset)
{
  if (is.null(mesh$normals))
    mesh <- vcgUpdateNormals(mesh)

  mesh$vb[1:3,] <- mesh$vb[1:3,]+t(offset*t(mesh$normals[1:3,]))
  invisible(mesh)
}
