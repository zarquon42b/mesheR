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

# export
oowriteh5.pPCA <- function(model,filename) {
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
    createH5Dataset(modinfo,"scores",t(model$PCA$x[,usePC]))
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
oh5model.read <- function(filename) {
    h5 <- H5File(filename)
    model <- list()
    class(model) <- "pPCA"
    mshape <- getH5Dataset(getH5Group(h5, "model"), "mean")[]
    model$mshape <- matrix(mshape,length(mshape)/3,3,byrow = T)
    modgroup <- getH5Group(h5, "model")
    Variance <- getH5Dataset(getH5Group(h5, "model"), "pcaVariance")[]
    model$PCA <- list()
    model$PCA$rotation <- getH5Dataset(getH5Group(h5, "model"), "pcaBasis")[]
    model$PCA$sdev <- sqrt(Variance)
    model$sigma <- getH5Dataset(modgroup, "noiseVariance")[]
    model <- setMod(model,model$sigma,exVar=1)
    modelinfogroup <- getH5Group(h5, "modelinfo")
    model$PCA$x <- getH5Dataset(modelinfogroup,"scores")
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

checkrhdf5 <- function() {
   # IP <- as.data.frame(installed.packages(),stringsAsFactors = F)$Package
    
    #check <- "rhdf5" %in% IP
    check <- require(rhdf5,quietly=T)
    if (!check) {
        ask <- readline("do you want to install missing package from Biconductor (internet required)? (y/n)\n")
    if (ask == "y") {
        source("http://bioconductor.org/biocLite.R")
        biocLite("rhdf5")
        require(rhdf5)
    } else
        stop("you will need to install rhdf5 first (Bioconductor only)")
    } else {
        require(rhdf5,quietly = T)
    }
 }
#' read statistical model as defined by statismo
#'
#' read statistical model as defined by statismo
#' @export
h5model.read <- function(filename) {
    checkrhdf5()
    model <- list()
    class(model) <- "pPCA"
    ##read model
    modread <- h5read(filename,"model",read.attributes = T)
    mshape <- modread$mean
    Variance <- as.vector(modread$pcaVariance)
    model$PCA <- list()
    model$PCA$rotation <- t(as.matrix(modread$pcaBasis))
    model$PCA$sdev <- sqrt(Variance)
    model$sigma <-  as.vector(modread$noiseVariance)
    ## read  modelinfo
    modelinfo <-  h5read(filename, "modelinfo")
    model$PCA$x <- as.matrix(modelinfo$scores)
    model$scale <- TRUE
    ## get representer
    representergroup <-  h5read(filename, "representer",read.attributes = T)
    referencetype <- as.character(attributes(representergroup)$datasetType)
    if (referencetype == "POLYGON_MESH"){
        refmesh <- list(); class(refmesh) <- "mesh3d"
        refmesh$vb <- rbind(t(representergroup$points),1)
        m <- ncol(refmesh$vb)-1
        refmesh$it <- t(representergroup$cells)+1
        if (nrow(refmesh$it) == 2)
            refmesh$it <- rbind(refmesh$it,refmesh$it[2,])
        model$refmesh <- refmesh
        k <- dim(refmesh$vb)
    } else {## to be filled with alternatives
        
    }
    model$mshape <- matrix(mshape,length(mshape)/3,3,byrow = T)
    model <- setMod(model,model$sigma,exVar=1)
    gc()
    H5garbage_collect()
    return(model)
    
}
#' @export
writeh5.pPCA <- function(model, filename) {
    checkrhdf5()
    if (file.exists(filename))
        file.remove(filename)
    h5createFile(filename)
    modgroup <- list()
    h5write(modgroup,filename,"model")
    mean <- as.vector(t(model$mshape))
    usePC <- model$usePC
    pcaBasis <- t(model$PCA$rotation[,usePC])
    pcaVariance <- as.vector(model$Variance[usePC,1])
    noiseVariance <- model$sigma
    h5createDataset(filename,"model/mean",dims=length(mean),H5type = "H5T_NATIVE_FLOAT")
    h5write(mean,filename,"model/mean")

    h5createDataset(filename,"model/pcaBasis",dims=dim(pcaBasis),H5type = "H5T_NATIVE_FLOAT")
    h5write(pcaBasis,filename,"model/pcaBasis")
    
    h5createDataset(filename,"model/noiseVariance",dims=1,H5type = "H5T_NATIVE_FLOAT")
    h5write(noiseVariance,filename,"model/noiseVariance")
    h5createDataset(filename,"model/pcaVariance",dims=length(pcaVariance),H5type = "H5T_NATIVE_FLOAT")
    h5write(pcaVariance,filename,"model/pcaVariance")
    
    modinfo <- list()
    h5write(modinfo,filename,"modelinfo")
    h5createDataset(filename,"modelinfo/scores",dims=dim(model$PCA$x),H5type = "H5T_NATIVE_FLOAT")
      
    h5write(model$PCA$x,filename,"modelinfo/scores")

    representer <- list()
    
    if (is.null(model$refmesh)) {
        attributes(representer) <- list(datasetType="POINT_SET",name="pointRepresenter")
        h5write(representer,filename,"representer")
        points <- model$mshape
        h5createDataset(filename,"representer/points",dims=dim(points),H5type = "H5T_NATIVE_FLOAT")
        h5write(points,filename,"representer/points")
        fid <- H5Fopen(filename)
        gid <- H5Gopen(fid,"representer")
        H5Gcreate(gid,"pointData")
        H5Gclose(gid)
        H5close()
        
    } else {
        attributes(representer) <- list(datasetType="POLYGON_MESH",name="meshRepresenter")
        h5write(representer,filename,"representer",write.attributes=T)
        points <- t(model$refmesh$vb[1:3,])
        cells <- t(model$refmesh$it-1)
        storage.mode(cells) <- "integer"
        h5createDataset(filename,"representer/points",dims=dim(points),H5type = "H5T_NATIVE_FLOAT")
        h5write(points,filename,"representer/points")
        h5createDataset(filename,"representer/cells",dims=dim(cells),storage.mode="integer")
        h5write(cells,filename,"representer/cells")
        fid <- H5Fopen(filename)
        gid <- H5Gopen(fid,"representer")
        H5Gcreate(gid,"pointData")
        H5Gcreate(gid,"cellData")
        H5Gclose(gid)
        H5close()
    }

    }






