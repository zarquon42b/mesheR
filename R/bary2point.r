#' transform barycentric coordinates into 3D-coordinates
#'
#' transform barycentric coordinates into 3D-coordinates
#' @param bary 3xn Matrix containing barycentric coordinates
#' @param faceptr integer vector of length n. Assigning a face to each triple of barycentric coordinates
#' @param mesh mesh to extract coordinates from
#'
#' @return a nx3 matrix containing 3D-coordinates
#' @examples
#' require(rgl)
#' require(Rvcg)
#' data(humface)
#' #extract 300 random points from humface 
#' coords <- vcgSample(humface,200)
#' #get barycentric coordinates
#' proj <- vcgClost(coords, humface, barycentric=TRUE)
#' #move original mesh
#' transface <- translate3d(humface, 10, 10 ,10)
#' newcoord <- bary2point(proj$barycoords, proj$faceptr, transface)
#' wire3d(transface, col=3)
#' spheres3d(newcoord,radius=0.5)
#' 
#' @export bary2point
bary2point <- function(bary,faceptr, mesh)
 {
     nbary <- ncol(bary)
     C <- Matrix::Matrix(0,nbary,ncol(mesh$vb))
     vertptr <- t(mesh$it[,faceptr])
     faceptr <- cbind(1:nbary,vertptr)
     for(i in 2:4)
         C[faceptr[,c(1,i)]] <- bary[i-1,]
     
     out <- as.matrix(C%*%t(mesh$vb[1:3,]))
     return(out)
 }

#' transfer points between two registered meshes
#'
#' transfer points between two registered meshes
#' @param x a matrix with 3D-coordinates (or a mesh), positioned on \code{mesh1}
#' @param mesh1 a mesh on which \code{x} is placed
#' @param mesh2 a mesh with vertices and faces corresponding to  \code{mesh1} (e.g. registered using \code{\link{gaussMatch}}).
#' @param tolwarn numeric: if at least one coordinate of \code{x} is further away from \code{mesh1} than \code{tolwarn}, a warning will be issued.
#' @details the function gets the barycentric coordinates of  \code{x} on \code{mesh1} and uses them to find the corresponding positions on \code{mesh2}
#' @return returns a matrix containing \code{x} tranfered to \code{mesh2}.
#' @examples
#' require(rgl)
#' require(Rvcg)
#' data(humface)
#' #extract 300 random points from humface 
#' coords <- vcgSample(humface,200,"pd")
#' 
#' #move original mesh
#' transface <- translate3d(humface, 10, 10 ,10)
#' ##extract coordinates
#' newcoord <- transferPoints(coords, humface, transface)
#' wire3d(transface, col=3)
#' spheres3d(newcoord)
#' @export
transferPoints <- function(x, mesh1, mesh2, tolwarn= 0.01) {
    
    chkv <- ncol(mesh1$vb) == ncol(mesh2$vb)
    chknf <- ncol(mesh1$it) == ncol(mesh2$it)
    if (!chkv)
        stop("amount of vertices does not correspond between meshes")
    if (!chknf)
        stop("amount of faces does not correspond between meshes")
    
    chkf <- which(mesh1$it != mesh2$it)
    if (length(chkf))
        stop("faces do not correspond between meshes")
    bary <- vcgClostKD(x,mesh1,k=10,barycentric=TRUE,sign=F)
    if (max(bary$quality) > tolwarn)
        warning(paste0("there are points that are farther away than ",tolwarn," mm"))
    out <- bary2point(bary$barycoord,bary$faceptr,mesh2)
    return(out)
}
        
     
    
