#' mix vertex colors of two registered meshes.
#' 
#' mixColorMesh mixes the vertex colors of two
#' registered meshes (identical faces and corresponding vertices). 
#' 
#' @param mesh1 triangular mesh of class "mesh3d".
#' @param mesh2 triangular mesh of class "mesh3d". 
#' @param alpha numeric: 0 <= alpha <=1 weight of color associated with mesh1.
#' @return returns mesh1 with mixed vertex colors.
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
#' @examples
#' require(rgl)
#' require(Morpho)
#' data(nose)
#' redmesh <- shortnose.mesh
#' redmesh$material$color <- matrix("#FF0000",dim(shortnose.mesh$it))
#' bluemesh <- shortnose.mesh
#' bluemesh$material$color <- matrix("#0000FF",dim(shortnose.mesh$it))
#' mixmesh <- mixColorMesh(bluemesh,redmesh)
#' shade3d(mixmesh)
#' @export mixColorMesh
mixColorMesh <- function(mesh1,mesh2,alpha=0.5)
    {
        col1 <- mesh1$material$color
        col2 <- mesh2$material$color

        nas1 <- which(is.na(col1))
        nas2 <- which(is.na(col2))
        if (length(nas1) > 0)
            col1[nas1] <- "#FFFFFF"
        if (length(nas2) > 0)
            col1[nas2] <- "#FFFFFF"

        mx <- hex(mixcolor(alpha=alpha,hex2RGB(col1),hex2RGB(col2)))
        mesh1$material$color <- matrix(mx,dim(mesh1$material$color))
        return(mesh1)
    }

