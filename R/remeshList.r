#' Remesh a list of registered meshes (with corresponding vertices)
#'
#' Remesh a list of registered meshes (with corresponding vertices), preserving correspondences using Quadric Edge Collapse decimation
#' @param matchlist a list of meshes with corresponding vertices (same amount and pseudo-homologous positions), e.g. registered with \code{\link{gaussMatch}}. Meshes do not need to be aligned.
#' @param reference integer: select the specimen the to which the decimation is to applied initially.
#' @param random if TRUE, a random specimen is selected for initial decimation
#' @param voxelSize voxel size for space discretization
#' @param discretize If TRUE, the position of the intersected edge of the marching cube grid is not computed by linear interpolation, but it is placed in fixed middle position. As a consequence the resampled object will look severely aliased by a stairstep appearance.
#' @param multiSample If TRUE, the distance field is more accurately compute by multisampling the volume (7 sample for each voxel). Much slower but less artifacts.
#' @param absDist If TRUE, an unsigned distance field is computed. In this case you have to choose a not zero Offset and a double surface is built around the original surface, inside and outside.
#' @param mergeClost logical: merge close vertices
#' @param silent logical: suppress messages
#' 
#' @return a list of remeshed meshes with correspondences preserved.
#' @details The remeshing is applied to a reference and then the barycentric coordinates of the new vertices on the original surface are calculated. These are used to extract the corresponding positions of the remeshed versions on all meshes in the sample. The remeshing is performed by the function \code{vcgUniformRemesh} from the Rvcg-package.
#' @importFrom Rvcg vcgUniformRemesh
#' @export
remeshList <- function(matchlist,reference=1,random=FALSE,voxelSize = NULL, discretize = FALSE, multiSample = FALSE, absDist = FALSE, mergeClost = FALSE,  silent = FALSE) {
    if (! random)
        reference <- sample(length(matchlist),size=1)

    ref <- matchlist[[reference]]
    remref <- vcgUniformRemesh(ref,voxelSize = voxelSize, discretize = discretize, multiSample = multiSample, absDist = absDist, mergeClost = mergeClost,  silent = silent)
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
#' decimate a list of registered meshes (with corresponding vertices)
#'
#' decimate a list of registered meshes (with corresponding vertices), preserving correspondences using Quadric Edge Collapse decimation
#'
#' @param matchlist a list of meshes with corresponding vertices (same amount and pseudo-homologous positions), e.g. registered with \code{\link{gaussMatch}}. Meshes do not need to be aligned.
#' @param reference integer: select the specimen the to which the decimation is to applied initially.
#' @param random if TRUE, a random specimen is selected for initial decimation
#' @param tarface Integer: set number of target faces.
#' @param percent Numeric: between 0 and 1. Set amount of reduction relative to
#' existing face number. Overrides tarface argument. 
#' @param edgeLength Numeric: tries to decimate according to a target mean edge
#' length. Under the assumption of regular triangles, the edges are half as
#' long by dividing the triangle into 4 regular smaller triangles.
#' @param topo logical: if TRUE, mesh topology is preserved.
#' @param quality logical: if TRUE, vertex quality is considered.
#' @param bound logical: if TRUE, mesh boundary is preserved.
#' @param optiplace logical: if TRUE, mesh boundary is preserved.
#' @param scaleindi logical: if TRUE, decimatiion is scale independent.
#' @param normcheck logical: if TRUE, normal directions are considered.
#' @param safeheap logical: if TRUE, safeheap update option enabled.
#' @param qthresh numeric: Quality threshold for decimation process.
#' @param boundweight numeric: Weight assigned to mesh boundaries.
#' @param normalthr numeric: threshold for normal check in radians.
#' @param silent logical, if TRUE no console output is issued.
#'
#' @return a list of decimated meshes with correspondences preserved.
#' @details The decimation is applied to a reference and then the barycentric coordinates of the new vertices on the original surface are calculated. These are used to extract the corresponding positions of the remeshed versions on all meshes in the sample. The Decimation is performed by the function \code{vcgQEdecim} from the Rvcg-package.
#' @importFrom Rvcg vcgQEdecim
#' @export
decimateList <- function(matchlist,reference=1,random=FALSE, tarface = NULL, percent = NULL, edgeLength = NULL, topo = FALSE, quality = TRUE, bound = FALSE, optiplace = TRUE,     scaleindi = TRUE, normcheck = FALSE, safeheap = FALSE, qthresh = 0.3,  boundweight = 1, normalthr = pi/2, silent = FALSE)  {
    if (!is.null(random))
        reference <- sample(length(matchlist),size=1)

    ref <- matchlist[[reference]]
    remref <- vcgQEdecim(ref,tarface = tarface, percent = percent, edgeLength = edgeLength, topo = topo, quality = quality, bound = bound, optiplace = optiplace,     scaleindi = scaleindi, normcheck = normcheck, safeheap = safeheap, qthresh = qthresh,  boundweight = boundweight, normalthr = boundweight, silent = silent)
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
