#' calculate or modify a conditional PCA based on 3D-coordinates
#'
#' calculate or modify a conditional PCA based on 3D-coordinates
#' 
#' @param array array of dimensions k x 3 x n, where k=number of coordinates and n=sample size.
#' @param missingIndex integer vector: specifies which points are missing in the conditional model
#' @param procMod object of class "procMod" as returned by pPCAcond or setMod
#' @param sigma estimate error variance (sensible is a value estimating coordinate error in terms of observer error)
#' @param exVar numeric value with \code{0 < exVar <= 1} specifying the PCs to be included by their cumulative explained Variance
#' @param refmesh a triangular mesh, where the vertices correspond to the coordinates in \code{array}
#' @param scale logical: allow scaling in Procrustes fitting
#' 
#' @return return a probabilistic PCA model of class "procMod"
#' 
#'
#' @examples
#' require(Morpho)
#' data(boneData)
#' model <- pPCAcond(boneLM[,,-1],missingIndex=3:4)
#' ## change parameters without recomputing Procrustes fit
#' model1 <- setMod(model, sigma=1, exVar=0.8)
#' 
#'
#'
#' @importFrom Morpho ProcGPA rotonmat arrMean3
#' @importFrom parallel mclapply detectCores
#' @rdname pPCAcond
#' @export
pPCAcond <- function(array, missingIndex, sigma=NULL, exVar=1,refmesh=NULL,scale=TRUE) {
    k <- dim(array)[1]
    use.lm=c(1:k)[-missingIndex]
    procMod <- ProcGPA(array[use.lm,,],scale=scale,CSinit=F,reflection=F,silent = TRUE)
    tmp <- array
    a.list <-  1:(dim(array)[3])
    tmp <- lapply(a.list, function(i) {mat <- rotonmat(array[,,i],array[use.lm,,i],procMod$rotated[,,i],scale=scale);return(mat)})
    tmp1 <- array
    for (i in 1:length(a.list))
        tmp1[,,i] <- tmp[[i]]

    procMod$rotated <- tmp1
    procMod$scale <- scale
    procMod$mshape <- arrMean3(tmp1)
    procMod$missingIndex <- missingIndex
    PCA <- prcomp(vecx(procMod$rotated,byrow = T))
    sel <- missingIndex*3
    sel <- sort(c(sel, sel-1, sel-2))
    sds <- PCA$sdev^2
    good <- which(sds > 1e-19)
    sds <- sds[good]
    PCA$rotation <- PCA$rotation[,good]
    PCA$sdev <- PCA$sdev[good]
    procMod$PCA <- PCA
    procMod$sel <- sel
    class(procMod) <- "pPCAcond"
    procMod <- setMod(procMod,sigma=sigma,exVar=exVar)
    procMod$refmesh <- refmesh
    return(procMod)
}
#' @rdname pPCAcond
#' @export
setMod <- function(procMod, sigma, exVar)UseMethod("setMod")

#' @rdname pPCAcond
#' @export
setMod.pPCAcond <- function(procMod,sigma=NULL,exVar=1) {
    k <- dim(procMod$rotated)[1]
    PCA <- procMod$PCA
    sds <- PCA$sdev^2
    sel <- procMod$sel
    sdsum <- sum(sds)
    sdVar <- sds/sdsum
    sdCum <- cumsum(sdVar)
    usePC <- which(sdCum <= exVar)
    Variance <- data.frame(eigenvalues=sds,exVar=sdVar, cumVar=sdCum)
    procMod$Variance <- Variance
    procMod$exVar <- exVar
    if (is.null(sigma))
        sigma <- 1/(3*k)*sum(sds[-usePC]) ##estimate sigma from remaining Variance

    if (sigma == 0)
        siginv <- 1e13
    else
        siginv <- sigma^-2

    sigest <- (sds - sigma^2)
    sigest <- sigest[which(sigest > 0)]
    usePC <- 1:min(length(usePC),length(sigest))
    procMod$usePC <- usePC
    W <- t(t(PCA$rotation[,usePC])*sqrt(sigest[usePC]))
    Wb <- W[-sel,]
    WbtWb <- crossprod(Wb)
    M <- siginv*WbtWb
    diag(M) <- diag(M)+1
    procMod$W <- W
    procMod$Wb <- Wb
    procMod$WbtWb <- WbtWb
    procMod$M <- M
    procMod$sigma <- sigma
    procMod$alphamean <- siginv*solve(M)%*%t(Wb)
    print(procMod,Variance=FALSE)
    return(procMod)
}

#' @export
print.pPCAcond <- function(x, digits = getOption("digits"), Variance=TRUE){
    cat(paste("   sigma =",x$sigma,"\n"))
    cat(paste("   exVar =",x$exVar,"\n\n"))
    cat(paste(" first",length(x$usePC),"of",ncol(x$PCA$rotation),"PCs used\n"))
    if (Variance) {
        cat("\n\n Model Variance:\n")
        print(x$Variance)
    }
}
#' @export
print.pPCA <- function(x, digits = getOption("digits"), Variance=TRUE){
    cat(paste("   sigma =",x$sigma,"\n"))
    cat(paste("   exVar =",x$exVar,"\n\n"))
    cat(paste(" first",length(x$usePC),"of",ncol(x$PCA$rotation),"PCs used\n"))
    if (Variance) {
        cat("\n\n Model Variance:\n")
        print(x$Variance)
    }
}
#' @export
restrictMissing <- function(x, model, refmesh,sdmax) UseMethod("restrictMissing")

restrictMissing.matrix <- function(x, model,refmesh=FALSE,sdmax) {
    
    mshape <- model$mshape
    missingIndex <- model$missingIndex
    rotsb <- rotonto(mshape[-missingIndex,],x,scale=model$scale)
    sb <- rotsb$yrot
    sbres <- sb-mshape[-missingIndex,]
    alpha <- model$alphamean%*%as.vector(t(sbres))
    if (!missing(sdmax)) {
        signalpha <- sign(alpha)
        alpha <- abs(alpha)
        outlier <- which(alpha > sdmax)
        alpha[outlier] <- sdmax
        alpha <- alpha*signalpha
    }
    
  
    
    #as.vector(W[,good]%*%alpha)
    estim <- t(as.vector(model$W%*%alpha)+t(mshape))
    estim <- rotreverse(estim,rotsb)
    if (!is.null(model$refmesh) && refmesh) {
        estimmesh <- model$refmesh
        estimmesh$vb[1:3,] <- t(estim)
        estimmesh <- vcgUpdateNormals(estimmesh)
        estim <- estimmesh
    }
    
    return(estim)
}
#restrictMissing.mesh3d <- function(x,model)
#    mat <- t(x$vb[1:3,])
#    estim <- restrictMissing(mat,model)
    

#' @export
pPCA <- function(array, sigma=NULL,exVar=1,scale=TRUE) {
    k <- dim(array)[1]
    procMod <- ProcGPA(array,scale=scale,CSinit=F,reflection=F)
    PCA <- prcomp(vecx(procMod$rotated,byrow = T))
    sds <- PCA$sdev^2
    good <- which(sds > 1e-19)
    sds <- sds[good]
    PCA$rotation <- PCA$rotation[,good]
    PCA$sdev <- PCA$sdev[good]
    procMod$PCA <- PCA
    procMod$scale <- scale
    class(procMod) <- "pPCA"
    procMod <- setMod(procMod,sigma=sigma,exVar=exVar)
    return(procMod)

}
#' @export
setMod.pPCA <- function(procMod,sigma=NULL,exVar=1) {
    k <- dim(procMod$rotated)[1]
    PCA <- procMod$PCA
    sds <- PCA$sdev^2
    sel <- procMod$sel
    sdsum <- sum(sds)
    sdVar <- sds/sdsum
    sdCum <- cumsum(sdVar)
    usePC <- which(sdCum <= exVar)
    Variance <- data.frame(eigenvalue=sds,exVar=sdVar, cumVar=sdCum)
    procMod$Variance <- Variance
    procMod$exVar <- exVar
    if (is.null(sigma))
        sigma <- 1/(3*k)*sum(sds[-usePC]) ##estimate sigma from remaining Variance

    if (sigma == 0)
        siginv <- 1e13
    else
        siginv <- sigma^-2

    sigest <- (sds - sigma^2)
    sigest <- sigest[which(sigest > 0)]
    usePC <- 1:min(length(usePC),length(sigest))
    procMod$usePC <- usePC
    procMod$sigma <- sigma
    W <- t(t(PCA$rotation[,usePC])*sqrt(sigest[usePC]))
    Win <- (t(PCA$rotation[,usePC])*(1/sqrt(sigest[usePC])))
    procMod$W <- W
    procMod$Win <- Win
    return(procMod)
}

#' @export
restrictpPCA <- function(x,model,sdmax=2) {
    mshape <- model$mshape
    rotsb <- rotonto(mshape,x,scale=model$scale)
    sb <- rotsb$yrot
    sbres <- sb-mshape
   # W <- model$W
    alpha <- model$Win%*%as.vector(t(sbres))
    signalpha <- sign(alpha)
    alpha <- abs(alpha)
    outlier <- which(alpha > sdmax)
    alpha[outlier] <- sdmax
    
    alpha <- alpha*signalpha
    
    res <- t(as.vector(model$W%*%alpha)+t(mshape))
    res <- rotreverse(res,rotsb)
    #res <- t(as.vector(W%*%t(W)%*%as.vector(t(sb)))+t(mshape))
    return(list(res=res,alpha=alpha))
}
    
