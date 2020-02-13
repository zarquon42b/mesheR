#' Transfers vertex colors between meshes.
#' 
#' colorTrans transfers the vertex color of mesh2 onto mesh1 by interpolating
#' vertex colors of closest face. 
#' 
#' @param mesh1 triangular mesh of class "mesh3d".
#' @param mesh2 triangular mesh of class "mesh3d". 
#' @param tol maximal distance of closest point included in color transfer.

#' @return returns the mesh1 with transfered vertex colors.
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' require(Morpho)
#' data(nose)
#' redmesh <- shortnose.mesh
#' #color mesh red
#' redmesh$material$color <- matrix("#FF0000",dim(shortnose.mesh$it))
#' # mesh without colors
#' nocolmesh <- shortnose.mesh
#' mixmesh <- colorTrans(nocolmesh, redmesh)
#' \dontrun{
#' require(rgl)
#' shade3d(mixmesh)
#' }
#' @importClassesFrom colorspace LAB RGB color
#' @importFrom colorspace RGB hex hex2RGB
#' @export colorTrans 
colorTrans <- function(mesh1, mesh2, tol=1)
    {
        if (is.null(mesh2$material$color))
            stop("no color to transfer")
        vn <- ncol(mesh2$vb)
        clost <- vcgClostKD(mesh1,mesh2,sign=F)
        col <-  rep("#FFFFFF", vn)
        tmp1 <- data.frame(it = as.vector(mesh2$it))
        tmp1$rgb <- as.vector(mesh2$material$color)
        tmp1 <- unique(tmp1)
        col[tmp1$it] <- tmp1$rgb
        colout <- hex2RGB(col,gamma=T)
        hits <- which(clost$quality <= tol)
        if (length(hits) == 0)
            stop("no hits within the given tolerance - increase tol!")
        colmesh1 <- matrix(0,ncol(mesh1$vb),4)
        colmesh1[hits,1] <- clost$faceptr[hits]
        colmesh1[hits,2:4] <- t(mesh2$it[,clost$faceptr[hits]])
        vert1col <- as(colout[colmesh1[hits,2],],"LAB")
        vert2col <- as(colout[colmesh1[hits,3],],"LAB")
        vert3col <- as(colout[colmesh1[hits,4],],"LAB")

        #vert1col <- as(vert1col,"LAB")
        #vert2col <- as(vert2col,"LAB")
        #vert3col <- as(vert3col,"LAB")
        outcol <- vert1col
        outcol@coords <- (vert1col@coords+vert2col@coords+vert3col@coords)/3
        #over <- which(outcol@coords > 1)
        #outcol@coords[over] <- 1

        goodcoll <- as(outcol,"RGB")
        blackit <- rep(1,ncol(mesh1$vb))
        colmat <- RGB(blackit,blackit,blackit)
        #colmat <- matrix(1,ncol(mesh1$vb),3)
        colmat@coords[hits,] <- goodcoll@coords
        #colmat <- rgb(colmat[, 1], colmat[, 2], colmat[,3], maxColorValue = max(colmat))
        over <- which(colmat@coords > 1)
        if (length(over) > 0)
            colmat@coords[over] <- 1
        colmat <- hex(colmat)
        
        material <- list()
        material$color <- colmat
        mesh1$material$color <- material$color
        #return(list(mesh=mesh1,clost=clost,colmat=colmat,outcol=outcol))
        return(mesh1)
    }
