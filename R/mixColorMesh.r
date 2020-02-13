#' mix vertex colors of two registered meshes.
#' 
#' mixColorMesh mixes the vertex colors of two
#' registered meshes (identical faces and corresponding vertices). 
#' 
#' @param mesh1 triangular mesh of class "mesh3d".
#' @param mesh2 triangular mesh of class "mesh3d". 
#' @param alpha numeric: 0 <= alpha <=1 weight of color associated with mesh2.
#' @return returns mesh1 with mixed vertex colors.
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' require(Morpho)
#' data(nose)
#' redmesh <- shortnose.mesh
#' redmesh$material$color <- rep("#FF0000",,ncol(shortnose.mesh$vb))
#' bluemesh <- shortnose.mesh
#' bluemesh$material$color <- rep("#0000FF",ncol(shortnose.mesh$vb))
#' mixmesh <- mixColorMesh(bluemesh,redmesh)
#' \dontrun{
#' require(rgl)
#' shade3d(mixmesh)
#' }
#' @importFrom colorspace mixcolor
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
        mesh1$material$color <- mx
        return(mesh1)
    }

