#' remove vertices within a convex surface mesh
#' 
#' remove vertices within a convex surface mesh by checking if they intersect
#' the mesh along their normals.
#' 
#' This only works for convex meshes, for concave substructures this might
#' result in unwanted removals.
#' 
#' @param mesh triangular mesh of class "mesh3d" 
#' @param mindist numeric: restrict to those vertices whose intersection point
#' is further away as mindist.
#' @param maxdist nummeric: restrict to those vertices whose intersection point
#' is closer away as maxdist.
#' @param explode logical: instead of surface normals use the direction outward
#' from a center of the mesh. Per default, the center of the meshes bounding
#' box is used. 
#' @param center numeric vector of length 3, used in the case "explode=TRUE",
#' if another center is required.
#' @return returns a cleaned mesh
#' @author Stefan Schlager
#' @seealso \code{\link{rmVertex}}
#' @keywords ~kwd1 ~kwd2
#' @export rmInternals
rmInternals <- function(mesh,mindist=0,maxdist=1e12,explode=FALSE,center=NULL) {

    
    if (explode)
      {
        if (is.null(center))
          center <- apply(meshcube(mesh),2,mean)
        mesh$normals[1:3,] <- mesh$vb[1:3,]-center
        mesh$normals[1:3,] <- apply(mesh$normals[1:3,],2,function(x){x <- x/(sqrt(sum(x^2)))})
      }
    meshoff <- meshOffset(mesh,0.005)
    check <- vcgRaySearch(meshoff,mesh)
    intern <- which(as.logical((check$quality * (check$distance > mindist)) * (check$quality * (check$distance < maxdist))))
     

    if (length(intern))
        out <- rmVertex(mesh,intern)
    else out <- mesh
    cat(paste("removed",length(intern),"vertices\n"))
    invisible(out)
  }

    
