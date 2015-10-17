#' get vertex indices of vertices outside a bounding box
#'
#' get vertex indices of vertices outside a bounding box
#'
#' @param x a k x 3 matrix or a triangular mesh of class "mesh3d"
#' @param corners a 8 x 3 matrix defining the corners of a bounding box, or a triangular mesh of class "mesh3d" representing the meshBox
#' @param outside logical: if TRUE the indices outside the box are returned, those inside
#' @return returns the indices of coordinates inside/outside the bounding box.
#' @examples
#' require(Rvcg)
#' require(Morpho)
#' 
#' data(humface)
#' bbox <- getMeshBox(humface,extend=-0.2)
#' bboxmesh <-  getMeshBox(humface,extend=-0.2,tri=TRUE)
#' \dontrun{
#' require(rgl)
#' spheres3d(bbox)
#' }
#' outside <- outsideBBox(humface,bbox)
#' humshrink <- rmVertex(humface,outside)
#' \dontrun{
#' wire3d(humshrink,col=3)#inside
#' wire3d(humface,col=2)#outside
#' shade3d(bboxmesh,alpha=0.4)
#' }
#' @export
outsideBBox <- function(x, corners, outside=TRUE) {
    if (!inherits(corners,"mesh3d"))
        box <- makeBox(corners)
    else
        box <- corners
    clost <- vcgClost(x,box,sign=T,facenormals=T)
    #x$normals[1:3,] <- x$vb[1:3,]-clost$vb[1:3,]
    #test <- normcheck(clost,x)
    test <- clost$quality
    if (outside)
        outside <- which(test < 0)
    else
        outside <- which(test >= 0)
    return(outside)
    
}

#' create a bounding box as mesh
#'
#' create a bounding box as mesh
#' @param corners a 8 x 3 matrix defining the corners of a bounding box.
#' @importFrom Morpho invertFaces
#' @importFrom Rvcg vcgRaySearch
#' @export
makeBox <- function(corners) {
    corners <- getMeshBox(corners)
    bbox4 <- list(vb = rbind(t(corners),1))
    it <- c(4,2,1,1,3,4,2,5,1,6,5,2,6,8,7,7,5,6,4,3,7,7,8,4,7,3,1,1,5,7,4,6,2,8,6,4)
    bbox4$it <- matrix(it,3,length(it)/3)
    class(bbox4) <- "mesh3d"
    bbox4 <- vcgUpdateNormals(bbox4)
    cbb <- cSize(vert2points(bbox4))
    cbboff <- cSize(vert2points(meshOffset(bbox4,1)))
    if (cbb > cbboff)
        bbox4 <- invertFaces(bbox4)
    
    return(bbox4)
}

#' compute the bounding box of a triangular mesh
#'
#' compute the bounding box of a triangular mesh
#' @param mesh triangular mesh of class "mesh3d"
#' @param tri logical: if TRUE, triangular mesh representing the bounding box is returned
#' @param extend numeric: offset extending the bbox (or shrinking if negative)
#' @param pca if TRUE, the box will be aligned by the principal axes.
#' @return returns a matrix with the corners of the box or a triangular mesh if tri=TRUE.
#' @examples
#' require(Rvcg)
#' data(humface)
#' bbox <- getMeshBox(humface)
#' \dontrun{
#' require(rgl)
#' spheres3d(bbox)
#' }
#' bbmesh <- getMeshBox(humface,tri=TRUE)
#' \dontrun{
#' shade3d(bbmesh,alpha=0.4)
#' wire3d(humface)
#' }
#' @export
getMeshBox <- function(mesh,tri=FALSE,extend=0,pca=TRUE) {
    if (inherits(mesh,"mesh3d"))
        vv2 <- vert2points(mesh)
    else
        vv2 <- mesh
    if (pca) {
        pp <- prcomp(vv2)
        ranges <- apply(pp$x,2,range)
    } else {
        ranges <- apply(vv2,2,range)
    }    
    ranges <- apply(ranges,2,extendrange,f=extend)
    mc <- as.matrix(expand.grid(ranges[,1],ranges[,2],ranges[,3]))
    if (pca)
        bbox <- sweep(mc%*%t(pp$rotation),2,-pp$center)
    else
        bbox <- mc
    if (tri)
        bbox <- makeBox(bbox)
    return(bbox)
}

#' remove parts of a mesh outside a bounding box defined by another mesh
#'
#' remove parts of a mesh outside a bounding box defined by another mesh
#' @param mesh1 reference mesh of which the bounding box will be computed
#' @param mesh2 mesh to be cropped
#' @param extend numeric: offset extending the bbox (or shrinking if negative)
#' @param pca if TRUE, the box will be aligned by the principal axes.
#' @param outside logical: if TRUE the part outside the box is removed
#' @return returns the cropped mesh2
#' @export 
cropOutsideBBox <- function(mesh1,mesh2,extend=0,pca=TRUE,outside=TRUE) {
    bbox <- getMeshBox(mesh1,tri=TRUE,extend=extend, pca=pca)
    outsidev <- outsideBBox(mesh2,bbox,outside=outside)
    if (length(outsidev))
        mesh2 <- Morpho::rmVertex(mesh2,outsidev)
    return(mesh2)
}
