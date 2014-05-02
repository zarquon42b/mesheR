#' calculate or modify a probablistic PCA based on 3D-coordinates
#'
#' calculate or modify a probablistic PCA based on 3D-coordinates
#' 
#' @encoding utf8
#' @param array array of dimensions k x 3 x n, where k=number of coordinates and n=sample size.
#' @param missingIndex integer vector: specifies which points are missing in the conditional model
#' @param procMod object of class "procMod" as returned by pPCAcond or setMod
#' @param sigma estimate error variance (sensible is a value estimating coordinate error in terms of observer error)
#' @param exVar numeric value with \code{0 < exVar <= 1} specifying the PCs to be included by their cumulative explained Variance
#' @param refmesh a triangular mesh, where the vertices correspond to the coordinates in \code{array}
#' @param scale logical: allow scaling in Procrustes fitting
#' @param x vector of deviation in standard deviations, coordinate matrix or triangular mesh of class "mesh3d" to be predicted.
#' @param sdmax maximum allowed standard deviation (per Principal axis) within the model space. Defines the probabilistic boundaries.
#' @param origSpace logical: rotate the estimation back into the original coordinate system.
#' @param pPCA logical: if TRUE, a constrained pPCA model is returned.
#' @param model probabilistic model of class "pPCA" or "pPCAcond"
#' @return \code{pPCA} and \code{pPCAcond} return a probabilistic PCA model of class "pPCA" or "pPCAcond" respectively. 
#' \code{predictPCA} and \code{predictPCAcond} select the most probable shape within a given model (within defined boundaries),
#' \code{setMod} is used to modify existing models without recomputing Procrustes registration and subsequent PCA.
#' @examples
#' require(Morpho)
#' data(boneData)
#' model <- pPCAcond(boneLM[,,-1],missingIndex=3:4)
#' ## change parameters without recomputing Procrustes fit
#' model1 <- setMod(model, sigma=1, exVar=0.8)
#' 
#'
#' @references
#' \enc{Lüthi}{Luethi} M, Albrecht T, Vetter T. 2009. Probabilistic modeling and visualization of the flexibility in morphable models. In: Mathematics of Surfaces XIII. Springer. p 251-264
#' 
#' @importFrom Morpho ProcGPA rotonmat arrMean3 vecx
#' @importFrom parallel mclapply detectCores
#' @importFrom Rvcg vcgUpdateNormals
#' @rdname pPCA
#' @export
pPCA <- function(array, sigma=NULL,exVar=1,scale=TRUE,refmesh=NULL) {
    k <- dim(array)[1]
    procMod <- ProcGPA(array,scale=scale,CSinit=F,reflection=F)
    PCA <- prcomp(vecx(procMod$rotated,byrow = T))
    sds <- PCA$sdev^2
    good <- which(sds > 1e-13)
    sds <- sds[good]
    PCA$rotation <- PCA$rotation[,good]
    PCA$sdev <- PCA$sdev[good]
    procMod$PCA <- PCA
    procMod$scale <- scale
    class(procMod) <- "pPCA"
    procMod$refmesh <- refmesh
    procMod <- setMod(procMod,sigma=sigma,exVar=exVar)
    return(procMod)

}
#' @rdname pPCA
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
    good <- which(sds > 1e-13)
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
#' @rdname pPCA
#' @export
setMod <- function(procMod, sigma, exVar)UseMethod("setMod")

#' #' @rdname pPCA
#' @export
setMod.pPCA <- function(procMod,sigma=NULL,exVar=1) {
    k <- dim(procMod$mshape)[1]
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

#' @rdname pPCA
#' @export
setMod.pPCAcond <- function(procMod,sigma=NULL,exVar=1) {
    k <- dim(procMod$mshape)[1]
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
    Minv <- solve(M)
    Minv <- (Minv+t(Minv))/2
    procMod$Minv <- solve(M)
    procMod$sigma <- sigma
    procMod$alphamean <- siginv*procMod$Minv%*%t(Wb)
    print(procMod,Variance=FALSE)
    return(procMod)
}

#' @export
print.pPCAcond <- function(x, digits = getOption("digits"), Variance=TRUE,...){
    cat(paste("   sigma =",x$sigma,"\n"))
    cat(paste("   exVar =",x$exVar,"\n\n"))
    cat(paste(" first",length(x$usePC),"of",ncol(x$PCA$rotation),"PCs used\n"))
    if (Variance) {
        cat("\n\n Model Variance:\n")
        print(x$Variance)
    }
}
#' @export
print.pPCA <- function(x, digits = getOption("digits"), Variance=TRUE,...){
    cat(paste("   sigma =",x$sigma,"\n"))
    cat(paste("   exVar =",x$exVar,"\n\n"))
    cat(paste(" first",length(x$usePC),"of",ncol(x$PCA$rotation),"PCs used\n"))
    if (Variance) {
        cat("\n\n Model Variance:\n")
        print(x$Variance)
    }
}
#' @rdname pPCA
#' @export
predictPCAcond <- function(x, model, refmesh,sdmax,pPCA=FALSE,...) UseMethod("predictPCAcond")

#' #' @rdname pPCA
#' @export
predictPCAcond.matrix <- function(x, model,refmesh=FALSE,sdmax,origSpace=TRUE,pPCA=FALSE,...) {
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
    ##as.vector(W[,good]%*%alpha)
    estim <- t(as.vector(model$W%*%alpha)+t(mshape))
    if (pPCA)
        procMod <- cond2pPCA(model,estim)

    if (origSpace)
        estim <- rotreverse(estim,rotsb)
    
    if (!is.null(model$refmesh) && class(model$refmesh) == "mesh3d" && refmesh) {
        estimmesh <- model$refmesh
        estimmesh$vb[1:3,] <- t(estim)
        estimmesh <- vcgUpdateNormals(estimmesh)
        estim <- estimmesh
    }
    if (pPCA)
        return(list(estim=estim,pPCA=procMod))
    else
        return(estim)
}

#' @rdname pPCA
#' @export
predictPCAcond.mesh3d <- function(x,model,refmesh=FALSE,sdmax, origSpace=TRUE,pPCA=FALSE,...) {
    mat <- t(x$vb[1:3,])
    estim <- predictPCAcond(x=mat,model=model,refmesh=refmesh,sdmax=sdmax,origSpace=origSpace)
    return(estim)
}


#' @rdname pPCA
#' @export
predictpPCA <- function(x,model,refmesh=FALSE,...)UseMethod("predictpPCA")

#' @rdname pPCA
#' @export
predictpPCA.matrix <- function(x,model,refmesh=FALSE,sdmax=2,origSpace=TRUE,...) {
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
    
    estim <- t(as.vector(model$W%*%alpha)+t(mshape))
    if (origSpace)
        estim <- rotreverse(estim,rotsb)

    if (!is.null(model$refmesh) && class(model$refmesh) == "mesh3d" && refmesh) {
        estimmesh <- model$refmesh
        estimmesh$vb[1:3,] <- t(estim)
        estimmesh <- vcgUpdateNormals(estimmesh)
        estim <- estimmesh
    }
    return(estim)
}

#' @rdname pPCA
#' @export
predictpPCA.mesh3d <- function(x,model,refmesh=FALSE,sdmax=2,origSpace=TRUE,...) {
    mat <- t(x$vb[1:3,])
    estim <- predictPCA(x=mat,model=model,refmesh=refmesh,sdmax=sdmax,origSpace=origSpace)
    return(estim)
}

#' @rdname pPCA
#' @export
predictpPCA.numeric <- function(x,model,refmesh=FALSE,...) {
    W <- model$W
    useit <- 1:(min(length(x),ncol(W)))
    if (length(useit) > 1)
        estim <- as.vector(W[,useit]%*%x)
    else
        estim <- as.vector(W[,useit]*x)
    
    estim <- t(estim+t(model$mshape))
    if (!is.null(model$refmesh) && class(model$refmesh) == "mesh3d" && refmesh) {
        estimmesh <- model$refmesh
        estimmesh$vb[1:3,] <- t(estim)
        estimmesh <- vcgUpdateNormals(estimmesh)
        estim <- estimmesh
    }
    return(estim)
}
    

cond2pPCA <- function(cpPCA, newMean) {
    procMod <- cpPCA
    procMod$mshape <- newMean
    eigM <- eigen(procMod$Minv, symmetric = T)
    sdev <- Re(eigM$values)
    good <- which(sdev > 1e-13)
    sdev <- sdev[good]
    procMod$usePC <- good
    procMod$PCA$sdev <- sqrt(sdev)
    newW <- procMod$W%*%Re(eigM$vectors[,good])
    procMod$W <- t(t(newW)*sdev)
    procMod$PCA$rotation <- newW
    class(procMod) <- "pPCA"
    allNames <- names(procMod)
    rem <- which(allNames %in% c("M","Minv","Wb","WbtWb","missingIndex","sel","alphamean"))
    procMod[rem] <- NULL
    procMod <- setMod(procMod,sigma=0,exVar=1)
    return(procMod)
}
