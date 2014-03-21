#' @export
h5model.write <- function(model,filename) {
    mean <- as.vector(t(model$mshape))
    pcaBasis <- model$PCs
    pcaVariance <- model$Variance[,1]
    h5 <- H5File(filename,'w')
    modgroup <- createH5Group(h5,"model")
    createH5Dataset(modgroup,"mean",mean)
    createH5Dataset(modgroup,"pcaBasis",pcaBasis)
    createH5Dataset(modgroup,"pcaVariance",pcaVariance)
    
}
h5model.read <- function(filename) {
    h5 <- H5File(filename)
    model <- list()
    mshape <- getH5Dataset(getH5Group(h5, "model"), "mean")[]
    model$mshape <- matrix(mshape,length(mshape),3,byrow = T)
    model$Variance <- getH5Dataset(getH5Group(h5, "model"), "pcaVariance")[]
    model$PCs <- getH5Dataset(getH5Group(h5, "model"), "pcaBasis")[]
    return(model)
    
}
