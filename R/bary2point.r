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
        #print(dim(C))
        return(out)
    }
