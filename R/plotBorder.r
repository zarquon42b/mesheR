#' plot all border edges of a triangular mesh
#' 
#' plot all border edges of a triangular mesh
#' 
#' 
#' @param mesh triangular mesh of class "mesh3d".  
#' @param col color of lines.  
#' @param lwd line width
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
#' @examples
#' require(rgl)
#' data(nose)
#' plotBorder(shortnose.mesh)
#' wire3d(shortnose.mesh,col=3)
#' @export plotBorder
plotBorder <- function(mesh,col=2,lwd=2)
    {
        edges <- vcgGetEdge(mesh)
        bord <- which(edges$border == 1)
        bordmesh <- mesh
        bordmesh$it <- rbind(edges$vert1[bord],edges$vert2[bord],edges$vert1[bord])
        bordmesh$material$color <- NULL
        wire3d(bordmesh,col=col,lwd=lwd,lit=F)
    }
