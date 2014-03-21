#' @export
h5model.write <- function(model,filename) {
    mean <- model$mshape
    pcaBasis <- model$PCs
    pcaVariance <- model$Variance[,1]
    h5 <- H5File(filename,'w')
    modgroup <- createH5Group(h5,"model")
    createH5Dataset(modgroup,"mean",mean)
    createH5Dataset(modgroup,"pcaBasis",pcaBasis)
    createH5Dataset(modgroup,"pcaVariance",pcaVariance)
    
}
