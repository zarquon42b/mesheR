#' consistently remesh a list of registered meshes (with corresponding vertices)
#'
#' consistently remesh a list of registered meshes (with corresponding vertices), preserving correspondences using Quadric Edge Collapse decimation
#'
#' @param matchlist a list of meshes with corresponding vertices (same amount and pseudo-homologous positions), e.g. registered with \code{\link{gaussMatch}}. Meshes do not need to be aligned.
#' @param reference integer: select the specimen the to which the decimation is to applied initially.
#' @param random if TRUE, a random specimen is selected for initial decimation
#' @param fun specify function to operate on the reference
#' @param ... parameters passed to fun
#' @return a list of decimated meshes with correspondences preserved.
#' @details The decimation is applied to a reference and then the barycentric coordinates of the new vertices on the original surface are calculated. These are used to extract the corresponding positions of the remeshed versions on all meshes in the sample. The Decimation is performed by the function \code{vcgQEdecim} from the Rvcg-package.
#' @importFrom Rvcg vcgQEdecim
#' @export
remeshList <- function(matchlist,reference=1,random=FALSE, fun=vcgQEdecim,...)  {

    if ((random))
        reference <- sample(length(matchlist),size=1)

    ref <- matchlist[[reference]]
    remref <- fun(ref,...)
    bary <- vcgClostKD(remref,ref,barycentric=T)

    outlist <- list()
    for (i in 1:length(matchlist)) {
        tmp <- remref
        tmpvert <- bary2point(bary$barycoords,bary$faceptr,matchlist[[i]])
        tmp$vb[1:3,] <- t(tmpvert)
        tmp <- vcgUpdateNormals(tmp)
        outlist[[i]] <- tmp
    }
    names(outlist) <- names(matchlist)
    return(outlist)
}
