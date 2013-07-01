#' move the centroid of a mesh to a given coordinate
#'
#' move the centroid of a mesh to a given coordinate
#' @param x triangular mesh of class 'mesh3d'
#' @param center coordinate where to translate the mesh, if center=NULL the centroid will be moved to the origin
#'
#' @return returns the translated mesh
#' @examples
#' data(humface)
#' wire3d(humface, col=3)
#' humcenter <- centermesh3d(humface)
#' #view translated mesh
#' wire3d(humcenter, col=2)
#' @export centermesh3d
centermesh3d <- function(x,center=NULL)
  {

    bary <- apply(vert2points(x),2,mean)
    if (is.null(center))
     center <- c(0,0,0)
     
       
    x$vb[1:3,] <- x$vb[1:3,]- (bary-center)
    invisible(x)
  }
 
