#' calculate or modify a probablistic PCA based on 3D-coordinates
#'
#' calculate or modify a probablistic PCA based on 3D-coordinates
#' 
#' @encoding utf8
#' @param array array of dimensions k x 3 x n, where k=number of coordinates and n=sample size.
#' @param missingIndex integer vector: specifies which points are missing in the conditional model
#' @param deselect logical: if TRUE, missingIndex references the existing coordinates instead of the missing ones.
#' @param procMod object of class "procMod" as returned by pPCAcond or setMod
#' @param sigma estimate error variance (sensible is a value estimating coordinate error in terms of observer error)
#' @param exVar numeric value with \code{0 < exVar <= 1} specifying the PCs to be included by their cumulative explained Variance
#' @param refmesh a triangular mesh, where the vertices correspond to the coordinates in \code{array}
#' @param scale logical: allow scaling in Procrustes fitting
#' @param fullfit logical: if FALSE only the non-missing points will be used for registration.
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
    procMod <- ProcGPA(array,scale=scale,CSinit=F,reflection=F) ##register all data using Procrustes fitting
    PCA <- prcomp(vecx(procMod$rotated,byrow = T),tol = sqrt(.Machine$double.eps)) ## calculate PCA
    sds <- PCA$sdev^2
    good <- which(sds > 1e-13)
    sds <- sds[good] ## remove PCs with very little variability
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
pPCAcond <- function(array, missingIndex,deselect=FALSE,sigma=NULL, exVar=1,refmesh=NULL,scale=TRUE,fullfit=FALSE) {
    k <- dim(array)[1]
    if (deselect) {
        missingIndex <- c(1:k)[-missingIndex]
    }
    use.lm=c(1:k)[-missingIndex]
    if (!fullfit) {
        procMod <- ProcGPA(array[use.lm,,],scale=scale,CSinit=F,reflection=F,silent = TRUE)##register all data using Procrustes fitting based on the non-missing coordinates
        tmp <- array
        a.list <-  1:(dim(array)[3])
        tmp <- lapply(a.list, function(i) {mat <- rotonmat(array[,,i],array[use.lm,,i],procMod$rotated[,,i],scale=scale,reflection = F);return(mat)})
        tmp1 <- array
        for (i in 1:length(a.list))
            tmp1[,,i] <- tmp[[i]]
        procMod$rotated <- tmp1
        procMod$mshape <- arrMean3(tmp1)
    } else {
        procMod <- ProcGPA(array,scale=scale,CSinit = F,reflection = F,silent = T)
    }
    procMod$use.lm <- use.lm
    procMod$scale <- scale
   
    procMod$missingIndex <- missingIndex
    PCA <- prcomp(vecx(procMod$rotated,byrow = T),tol = sqrt(.Machine$double.eps)) ## calculate PCA
    sel <- missingIndex*3
    sel <- sort(c(sel, sel-1, sel-2))
    sds <- PCA$sdev^2
    good <- which(sds > 1e-13)
    sds <- sds[good]## remove PCs with very little variability
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

#' @rdname pPCA
#' @export
setMod.pPCA <- function(procMod,sigma=NULL,exVar=1) {
    k <- dim(procMod$mshape)[1]
    PCA <- procMod$PCA
    sds <- PCA$sdev^2
    sdsum <- sum(sds)
    sdVar <- sds/sdsum
    sdCum <- cumsum(sdVar)
    if (is.null(sigma))
        sigma <- 1/(3*k)*sum(sds[-usePC]) ##estimate sigma from remaining Variance

    if (sigma >= sdsum) {
        warning(paste("sigma > overall variance set to",sdsum/2))
        sigma <- sdsum/2
    }
    usePC <- which(sdCum <= exVar)
    
        
    Variance <- data.frame(eigenvalue=sds,exVar=sdVar, cumVar=sdCum) ##make Variance table 
    procMod$Variance <- Variance
    
    if (sigma == 0)
        siginv <- 1e13
    else
        siginv <- sigma^-2

    sigest <- (sds - sigma^2)
    sigest <- sigest[which(sigest > 0)]
    usePC <- 1:max(1,min(length(usePC),length(sigest)))
    procMod$usePC <- usePC
    procMod$exVar <- sdCum[max(usePC)]##calculate variance explained by that model compared to that of the training sample
    procMod$sigma <- sigma
    W <- t(t(PCA$rotation[,usePC])*sqrt(sigest[usePC])) ##Matrix to project scaled PC-scores back into the config space
    Win <- (t(PCA$rotation[,usePC])*(1/sqrt(sigest[usePC]))) ##Matrix to project from config space into the scaled PC-space
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

    if (is.null(sigma))
        sigma <- 1/(3*k)*sum(sds[-usePC]) ##estimate sigma from remaining Variance

    if (sigma >= sdsum) {
        warning(paste("sigma > overall variance set to",sdsum/2))
        sigma <- sdsum/2
    }
    usePC <- which(sdCum <= exVar)
    Variance <- data.frame(eigenvalues=sds,exVar=sdVar, cumVar=sdCum) ##make Variance table
    procMod$Variance <- Variance
    procMod$exVar <- exVar
   

    if (sigma == 0)
        siginv <- 1e13
    else
        siginv <- sigma^-2

    sigest <- (sds - sigma^2)
    sigest <- sigest[which(sigest > 0)]
    usePC <- 1:min(length(usePC),length(sigest))
    procMod$usePC <- usePC
    procMod$exVar <- sdCum[max(usePC)]

    W <- t(t(PCA$rotation[,usePC])*sqrt(sigest[usePC]))
    Wb <- W[-sel,]
    WbtWb <- crossprod(Wb)
    M <- siginv*WbtWb
    diag(M) <- diag(M)+1
    procMod$W <- W ##Matrix to project scaled PC-scores back into the config space
    procMod$Wb <- Wb
    procMod$WbtWb <- WbtWb
    procMod$M <- M
    stry <- try(Minv <- solve(M)) 
    if (inherits(stry,"try-error")) {
        Minv <- Morpho:::armaGinv(M)
        message("singular Matrix")
    }
    Minv <- (Minv+t(Minv))/2
    procMod$Minv <- Minv ## covariance structure of the alpha under the restriction based on non-missing data.
    procMod$sigma <- sigma
    procMod$alphamean <- siginv*procMod$Minv%*%t(Wb) ## the general mean of the constrained distribution
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
#' predict or restrict a mesh or matrix based on a statistical model
#'
#' predict or restrict a mesh or matrix based on a statistical model
#' @param x a matrix, a mesh3d or a vector (for pPCA models) containing standardized variables within the PC-space
#' @param model model of class pPCA or pPCAcond
#' @param refmesh if TRUE and the model contains a representer mesh, a surface mesh will be returned, coordinate matrix otherwise.
#' @param origSpace logical: rotate the estimation back into the original coordinate system.
#' @param pPCA logical: if TRUE, a constrained pPCA model is returned.
#' "chisq" uses the Chi-Square distribution of the squared Mahalanobisdistance, while "dist" restricts the values to
#' be within a multi-dimensional sphere of radius \code{sdmax}. If FALSE the probability will be determined per PC separately.
#' @param use.lm optional: integer vector specifying row indices of the coordinates to use for rigid registration on the model's meanshape.
#' @param sdmax maximum allowed standard deviation (per Principal axis) within the model space. Defines the probabilistic boundaries.
#' @param mahaprob character: if != "none", use mahalanobis-distance to determine overall probability (of the shape projected into the model space.
#' @return \code{predictpPCA} returns a matrix/mesh3d restricted to the boundaries given by the modelspace.
#'
#' \code{predictPCAcond} returns a list with
#' \item{estim}{matrix/mesh3d representing the mean of the restricted space}
#'
#' \item{pPCA}{if \code{pPCA = TRUE} a pPCA model representing the gaussian subspace given the constraints is returned}
#' \item{rot}{the transformation of x into the modelspace that can be reverted by calling \code{rotreverse} from the package Morpho} 

#' @rdname predictpPCA
#' @export
predictPCAcond <- function(x, model, refmesh, origSpace=TRUE, pPCA=FALSE,...) UseMethod("predictPCAcond")

#' @rdname predictpPCA
#' @export
predictPCAcond.matrix <- function(x, model,refmesh=FALSE,origSpace=TRUE,pPCA=FALSE,...) {
    mshape <- model$mshape
    missingIndex <- model$missingIndex
    use.lm <- model$use.lm
    rotsb <- rotonto(mshape[use.lm,],x,scale=model$scale,reflection = F)
    sb <- rotsb$yrot
    sbres <- sb-mshape[use.lm,]
    alpha <- model$alphamean%*%as.vector(t(sbres))
    
   
    ##as.vector(W[,good]%*%alpha)
    estim <- t(as.vector(model$W%*%alpha)+t(mshape))
    if (pPCA)
        procMod <- as.pPCA(model,estim)

    if (origSpace)
        estim <- rotreverse(estim,rotsb)
    
    if (!is.null(model$refmesh) && class(model$refmesh) == "mesh3d" && refmesh) {
        estimmesh <- model$refmesh
        estimmesh$vb[1:3,] <- t(estim)
        estimmesh <- vcgUpdateNormals(estimmesh)
        estim <- estimmesh
    }
    if (pPCA)
        return(list(estim=estim,pPCA=procMod,rot=rotsb))
    else
        return(estim)
}

#' @rdname predictpPCA
#' @export
predictPCAcond.mesh3d <- function(x,model,refmesh=FALSE, sdmax, origSpace=TRUE,pPCA=FALSE,...) {
    mat <- t(x$vb[1:3,])
    estim <- predictPCAcond(x=mat,model=model,refmesh=refmesh,sdmax=sdmax,origSpace=origSpace)
    return(estim)
}


#' @rdname predictpPCA
#' @export
predictpPCA <- function(x,model,refmesh=FALSE,...)UseMethod("predictpPCA")

#' @rdname predictpPCA
#' @export
predictpPCA.matrix <- function(x,model,refmesh=FALSE,origSpace=TRUE,use.lm=NULL,sdmax,mahaprob=c("none","chisq","dist"),...) {
    mahaprob <- substr(mahaprob[1],1L,1L)
    mshape <- model$mshape
    if (is.null(use.lm)) {
        rotsb <- rotonto(mshape,x,scale=model$scale,reflection = F)
        sb <- rotsb$yrot
    } else {
        rotsb <- rotonto(mshape[use.lm,],x[use.lm,],scale=model$scale,reflection=F)
        sb <- rotonmat(x,x[use.lm,],rotsb$yrot)
    }
    sbres <- sb-mshape
   # W <- model$W
    alpha <- model$Win%*%as.vector(t(sbres))
    sdl <- nrow(model$Win)
    
     if (!missing(sdmax)) {
         if (mahaprob != "n") {
            sdl <- nrow(model$W)
            probs <- sum(alpha^2)
            if (mahaprob == "c") {
                Mt <- qchisq(1-2*pnorm(sdmax,lower.tail=F),df=sdl)
                probs <- sum(alpha^2)
            } else if (mahaprob == "d") {
                Mt <- sdmax
                probs <- sqrt(probs)
            }
            if (probs > Mt ) {
                sca <- Mt/probs
                alpha <- alpha*sca
            }
        } else { 
            signalpha <- sign(alpha)
            alpha <- abs(alpha)
            outlier <- which(alpha > sdmax)
            alpha[outlier] <- sdmax
            alpha <- alpha*signalpha
        }
    }
    
    if ("coeffs" %in% names(list(...))) {
        return(alpha)
        
    } else {
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
}

#' @rdname predictpPCA
#' @export
predictpPCA.mesh3d <- function(x,model,refmesh=FALSE,origSpace=TRUE,use.lm=NULL,sdmax,mahaprob=c("none","chisq","dist"),...) {
    mat <- t(x$vb[1:3,])
    estim <- predictpPCA(vert2points(x),model=model,refmesh=refmesh,sdmax=sdmax,origSpace=origSpace,use.lm=use.lm,mahaprob=mahaprob,...)
    return(estim)
}

#' @rdname predictpPCA
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

#' @export
as.pPCA <- function(x,..)UseMethod("as.pPCA")

#' @export
as.pPCA.pPCAcond <- function(x, newMean,...) { #convert a pPCAcond to a pPCA by adding a new mean config # not to be called directly
    procMod <- x
    procMod$mshape <- newMean
    sds <- procMod$PCA$sdev^2
    usePC <- procMod$usePC
    sigma <- procMod$sigma
    sigest <- (sds[usePC] - sigma^2)
    Minv <- procMod$Minv
    udut <- t(t(Minv)*sigest)
    eigM <- eigen(udut,symmetric = T)
    sds <- Re(eigM$values)
    good <- which(sds > 1e-15)
    sds <- sds[good]
    procMod$usePC <- good
    procMod$PCA$sdev <- sqrt(sds)
    newW <- procMod$PCA$rotation[,usePC]%*%Re(eigM$vectors[,good])
    procMod$W <- t(t(newW)*sqrt(sds))
    procMod$Win <- t(newW)/sqrt(sds)
    procMod$PCA$rotation <- newW
    class(procMod) <- "pPCA"
    allNames <- names(procMod)
    rem <- which(allNames %in% c("M","Minv","Wb","WbtWb","missingIndex","sel","alphamean"))
    procMod[rem] <- NULL
    sdsum <- sum(sds)
    sdVar <- sds/sdsum
    sdCum <- cumsum(sdVar)
    Variance <- data.frame(eigenvalues=sds,exVar=sdVar, cumVar=sdCum)
    procMod$Variance <- Variance
    return(procMod)
}

#' @export
as.pPCAcond <- function(x, missingIndex,...)UseMethod("as.pPCAcond")#convert a pPCA to a pPCAcond by providing constrains # not to be called directly

#' @export
as.pPCAcond.pPCA <- function(x,missingIndex) {
    procMod <- x
    procMod$missingIndex <- missingIndex
    sel <- missingIndex*3
    sel <- sort(c(sel, sel-1, sel-2))
    procMod$sel <- sel
    class(procMod) <- "pPCAcond"
    procMod <- setMod(procMod,sigma=procMod$sigma,exVar=procMod$exVar)
    return(procMod)
}

#' calculate probability/coefficients for a matrix/mesh given a statistical model
#'
#' calculate probability for a matrix/mesh given a statistical model
#' @param x matrix or mesh3d
#' @param model a model of class pPCA
#' @param use.lm integer vector specifying row indices of the coordinates to use for rigid registration on the model's meanshape.
#' @return \code{getProb} returns a probability, while \code{getCoefficients} returns the (scaled) scores in the pPCA space.
#' @export
getProb <- function(x,model,use.lm) UseMethod("getProb")

#' @rdname getProb
#' @export
getProb.matrix <- function(x,model,use.lm=NULL) {
    mshape <- model$mshape
    if (is.null(use.lm)) {
        rotsb <- rotonto(mshape,x,scale=model$scale,reflection = F)
        sb <- rotsb$yrot
    } else {
        rotsb <- rotonto(mshape[use.lm,],x[use.lm,],scale=model$scale,reflection=F)
        sb <- rotonmat(x,x[use.lm,],rotsb$yrot)
    }
    sbres <- sb-mshape
   # W <- model$W
    alpha <- model$Win%*%as.vector(t(sbres))
    sdl <- nrow(model$Win)
    probs <- sum(alpha^2)
    probout <- pchisq(probs,lower.tail = F,df=sdl)
    return(probout)
}

#' @rdname getProb
#' @export
getProb.mesh3d <- function(x,model,use.lm=NULL) {
    x <- vert2points(x)
    out <- getProb(x,model=model,use.lm=use.lm)
    return(out)
}

#' @rdname getProb
#' @export
getCoefficients <- function(x, model,use.lm=NULL) {
    out <- predictpPCA(x,model,use.lm,coeffs=NULL)
    return(out)
}
    
