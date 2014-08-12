#' close two meshes where one is only an offset of the other
#'
#' close two meshes where one is only an offset of the other
#'
#' @param mesh1 mesh
#' @param mesh2 mesh
#' @param invert logical: if TRUE the orientation of mesh1 will be flipped
#' @importFrom Morpho mergeMeshes conv2backf
#' @importFrom Rvcg vcgGetEdge vcgClean
#' @export
closeMeshes <- function(mesh1,mesh2,invert=TRUE) {
    if (invert)
        mesh1 <- conv2backf(mesh1)
    nvb <- ncol(mesh1$vb)
    merge <- mergeMeshes(mesh1,mesh2)
    gE <- vcgGetEdge(mesh1)
    gE <- gE[gE$border == T,]
    itnew <- rbind(gE$vert1,gE$vert2,gE$vert1+nvb)
    itnew <- cbind(itnew,rbind(gE$vert1+nvb,gE$vert2+nvb,gE$vert2))
    merge$it <- cbind(merge$it,itnew)
    merge <- vcgClean(merge,sel=7)
    merge <- vcgUpdateNormals(merge)
    return(merge)
}
