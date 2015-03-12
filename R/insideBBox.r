#' get vertex indices of vertices outside a bounding box
#'
#' get vertex indices of vertices outside a bounding box
#'
#' @param x a k x 3 matrix or a triangular mesh of class "mesh3d"
#' @param corners a 8 x 3 matrix defining the corners of a bounding box.
#' @param outside logical: if TRUE the indices outside the box are returned, those insde
#' @return returns the indices of coordinates inside/outside the bounding box.
#' @examples
#' require(Rvcg)
#' require(Morpho)
#' require(rgl)
#' data(humface)
#' bbox <- getMeshBox(humface,extend=-0.2)
#' bboxmesh <-  getMeshBox(humface,extend=-0.2,tri=TRUE)
#' spheres3d(bbox)
#' outside <- outsideBBox(humface,bbox)
#' humshrink <- rmVertex(humface,outside)
#' wire3d(humshrink,col=3)
#' wire3d(humface,col=2)
#' shade3d(bboxmesh,alpha=0.2)
#' 
#' @export
outsideBBox <- function(x, corners, outside=TRUE) {
    box <- makeBox(corners)
    clost <- vcgClost(x,box,sign=T)
    if (outside)
        outside <- which(clost$quality > 0)
    else
        outside <- which(clost$quality < 0)
    return(outside)
    
}

#' create a bounding box as mesh
#'
#' create a bounding box as mesh
#' @param corners a 8 x 3 matrix defining the corners of a bounding box.
#' @export
makeBox <- function(corners) {
    bbox4 <- list(vb = rbind(t(corners),1))
    it <- c(1,2,4,4,3,1,1,5,2,2,5,6,7,8,6,6,5,7,7,3,4,4,8,7,1,3,7,7,5,1,2,6,4,4,6,8)
    
    bbox4$it <- matrix(it,3,length(it)/3)
    class(bbox4) <- "mesh3d"
    return(bbox4)
}

#' compute the bounding box of a triangular mesh
#'
#' compute the bounding box of a triangular mesh
#' @param mesh triangular mesh of class "mesh3d"
#' @param tri logical: if TRUE, triangular mesh representing the bounding box is returned
#' @param extend numeric: offset extending the bbox (or shrinking if negative)
#' @return returns a matrix with the corners of the box or a triangular mesh if tri=TRUE.
#' @examples
#' require(Rvcg)
#' require(rgl)
#' data(humface)
#' bbox <- getMeshBox(humface)
#' spheres3d(bbox)
#' bbmesh <- getMeshBox(humface,tri=TRUE)
#' wire3d(bbmesh)
#' wire3d(humface)
#' @export
getMeshBox <- function(mesh,tri=FALSE,extend=0) {
    vv2 <- vert2points(mesh)
    pp <- prcomp(vv2)
    ranges <- apply(pp$x,2,range)
    ranges <- apply(ranges,2,extendrange,f=extend)
    mc <- as.matrix(expand.grid(ranges[,1],ranges[,2],ranges[,3]))
    bbox <- sweep(mc%*%t(pp$rotation),2,-pp$center)
    if (tri)
        bbox <- makeBox(bbox)
    return(bbox)
}

