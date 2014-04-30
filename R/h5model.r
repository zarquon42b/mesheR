#' Export a probabilistic PCA model into the statismo format
#'
#' Export a probabilistic PCA model into the statismo format
#'
#' @param model PCA model
#' @param filename character specifying filename
#'
#' @references
#' https://docs.google.com/document/d/1ySwSsD-di8ez5WnViLl4VsYeRI-ggVa-HQRxG4kLjH4/edit
#' @export
writeh5 <- function(model,filename)UseMethod("writeh5")

#' @export
writeh5.pPCA <- function(model,filename) {
    if (file.exists(filename))
        file.remove(filename)
    mean <- as.vector(t(model$mshape))
    usePC <- model$usePC
    pcaBasis <- model$PCA$rotation[,usePC]
    pcaVariance <- model$Variance[usePC,1]
    noiseVariance <- model$sigma
    h5 <- H5File(filename,'w')
    modgroup <- createH5Group(h5,"model")
    createH5Dataset(modgroup,"mean",mean)
    createH5Dataset(modgroup,"pcaBasis",pcaBasis)
    createH5Dataset(modgroup,"pcaVariance",pcaVariance)
    createH5Dataset(modgroup,"noiseVariance",noiseVariance)
    modinfo <- createH5Group(h5,"modelinfo")
    createH5Dataset(modinfo,"scores",model$PCA$x[,usePC])
    representergroup <- createH5Group(h5,"representer")
    
    if (is.null(model$refmesh)) {
        createH5Attribute(representergroup,"name","pointRepresenter")
        pointgroup <- createH5Group(representergroup,"pointData")
        createH5Attribute(representergroup,"datasetType","POINT_SET")
        createH5Dataset(representergroup,"points",t(mean))
    } else {
        pointgroup <- createH5Group(representergroup,"pointData")
        cellgroup <- createH5Group(representergroup,"cellData")
        createH5Attribute(representergroup,"name","meshRepresenter")
        createH5Attribute(representergroup,"datasetType","POLYGON_MESH")
        createH5Dataset(representergroup,"points",model$refmesh$vb[1:3,])
        createH5Dataset(representergroup,"cells",model$refmesh$it-1)
    }
        
}
h5model.read <- function(filename) {
    h5 <- H5File(filename)
    model <- list()
    class(model) <- "pPCA"
    mshape <- getH5Dataset(getH5Group(h5, "model"), "mean")[]
    model$mshape <- matrix(mshape,length(mshape)/3,3,byrow = T)
    Variance <- getH5Dataset(getH5Group(h5, "model"), "pcaVariance")[]
    model$PCA <- list()
    model$PCA$rotation <- getH5Dataset(getH5Group(h5, "model"), "pcaBasis")[]
    model$PCA$sdev <- sqrt(Variance)
    model$sigma <- getH5Dataset(getH5Group(h5, "model"), "noiseVariance")[]
    model <- setMod(model,model$sigma,exVar=1)
    model$scale <- TRUE
    referencetype <- getH5Attribute(getH5Group(h5, "representer"), "datasetType")[]
    if (referencetype == "POLYGON_MESH"){
        refmesh <- list(); class(refmesh) <- "mesh3d"
        refmesh$vb <- rbind(getH5Dataset(getH5Group(h5, "representer"), "points")[],1)
        refmesh$it <- getH5Dataset(getH5Group(h5, "representer"), "cells")[]+1
        model$refmesh <- refmesh
    }        
    return(model)
    
}
