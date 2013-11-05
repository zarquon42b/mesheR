glView <- function()
  {
    out <- c(0,0,0)
    info <- rgl.projection()
    mdl <- info$model
    out[1] <- -(mdl[1] * mdl[13] + mdl[2] * mdl[14] + mdl[3] * mdl[15])
    out[2] <-  -(mdl[5] * mdl[13] + mdl[6] * mdl[14] + mdl[7] * mdl[15])
    out[3] <- -(mdl[9] * mdl[13] + mdl[10] * mdl[14] + mdl[11] * mdl[15])
    return(out)
  }


#' determine visibility of mesh's vertices
#' 
#' Determine which vertices of a rendered mesh are visible from the current
#' viewpoint
#' 
#' 
#' @param mesh triangular mesh of class "mesh3d". Must be currently rendered in
#' an rgl window.
#' @param offset initial offset to move vertex slightly away from the surface.
#' @return returns logical vector, assigning TRUE/FALSE to each vertex of a
#' mesh.
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
#' @seealso \code{\link{selectVertex}}, \code{\link{cutMesh}}
#' @examples
#' require(rgl)
#' require(Morpho)
#' data(nose)
#' shade3d(shortnose.mesh,col=3)
#' visi <- glVisible(shortnose.mesh)
#' points3d(vert2points(shortnose.mesh)[which(visi),])
#' @importFrom rgl rgl.projection
#' @export glVisible
glVisible <- function(mesh, offset=1e-3)
{
  mesh <- adnormals(mesh)
  mesh0 <- meshOffset(mesh, offset)
  viewpoint <- c(glView(),0)
  normals <- viewpoint-mesh0$vb
  #normals <- apply(normals,2,function(x) x <- x/sqrt(sum(x^2)))
  mesh0$normals <- normals
  tmp <- as.logical(vcgRaySearch(mesh0,mesh)$quality)
  out <- tmp
  out[tmp] <- FALSE
  out[!tmp] <- TRUE
  return(out)
}

    
