objectiveMSQ <- function(x,clost,A,B,tarclost) {
    out <- 0
    if (!is.null(A)){
        Ax <- A%*%x
        out <- sum((Ax-clost)^2)/(nrow(Ax)/3)
    }
    if (!is.null(B)) {
        Bx <- B%*%x
        out <- out+sum((Bx-tarclost)^2)/(nrow(Bx)/3)
    }
    return(out)
}


objectiveMSQ.grad <- function(x,clost,A,B,tarclost) {
    grad <- x*0
    if  (!is.null(A)){
        Ax <- A%*%x
        Axb <- Ax-clost
        grad <- (2*t(A)%*%Axb)/(nrow(Ax)/3)
    }
    if (!is.null(B)) {
        Bx <- B%*%x
        Bxb <- Bx-tarclost
        grad <- grad + (2*t(B)%*%Bxb)/(nrow(Bx)/3)
    }
    return(grad)
}

#' fit a model minimizing the (symmetric) mean squared  distance
#'
#' fit a model minimizing the (symmetric) mean squared distance
#' @param model statistical model of class pPCA
#' @param tarmesh a target mesh already aligned to the model
#' @param iterations numbers of iterations to run
#' @param lbfgs.iter integer:
#' @param lbfgs.orthantwise_c integer or integer vector of length iterations. Contains coefficient(x) for the \code{L1} norm of variables. This parameter should be set to zero for standard minimization
#' problems. Setting this parameter to a positive value
#' activates Orthant-Wise Limited-memory Quasi-Newton (OWL-QN)
#' method, which minimizes the objective function \code{F(x)}
#' combined with the L1 norm \code{|x|} of the variables, \code{{F(x) + C
#' |x|}}. This parameter is the coefficient for the \code{|x|}, i.e.,
#' \code{C}. As the L1 norm \code{|x|} is not differentiable at zero, the
#' library modifies function and gradient evaluations from a
#' client program suitably. The default value is zero. Note that
#' the objective function minimized by alternative packages
#' (e.g., \code{glmnet}) is of the form : \code{F(x)/N + C |x|}, where \code{N}
#' is the number of parameters. \code{lbfgs} does not divide the
#' likelihood function by \code{N}. To achieve equivalence with
#' \code{glmnet} result, take this difference of implementation into
#' account.
#' @param refdist maximal distance from model to reference to be considered
#' @param tardist maximal distance from target to model to be considered
#' @param rho numeric: allowed normal deviation of a point to be considered as corresponding.
#' @param symmetric integer: specify which MSE to minimize. 0=search both ways, 1=model to target, 2=target to model.
#' @param sdmax constrain parameters (normalized PC-scores) to be within +- sdmax
#' @param mahaprob character: if != "none", use mahalanobis-distance to determine overall probability (of the shape projected into the model space."chisq" uses the Chi-Square distribution of the squared Mahalanobisdistance, while "dist" restricts the values to be within a multi-dimensional sphere of radius \code{sdmax}. If FALSE the probability will be determined per PC separately.
#' @param initparams a vector with initial estimations of the model parameters. Set to zeros if NULL.
#' @param k integer: amount of closest faces to consider during closest point search.
#' @param threads integer: number of cores to use for closest point search
#' @param method optimizer method. Can be one of "lbfgs", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent". lbfgs calls the function from package lbfgs, while the others call \code{\link{optim}}.
#' @param silent logical: if TRUE output will be suppressed.
#' @param ... additional parameters to be passed to \code{\link{lbfgs}} or \code{\link{optim}}.
#' @return
#' \item{par}{the model's parameters}
#' \item{mesh}{the fitted mesh}
#' @examples
#' \dontrun{
#' require(RvtkStatismo)
#' download.file(url="https://github.com/marcelluethi/
#' statismo-shaperegistration/raw/master/data/VSD001_femur.vtk","./VSD001_femur.vtk",method = "w")
#' download.file(url="https://github.com/marcelluethi/
#' statismo-shaperegistration/raw/master/data/VSD002_femur.vtk","./VSD002_femur.vtk",method = "w")
#' download.file(url="https://github.com/marcelluethi/
#' statismo-shaperegistration/raw/master/data/VSD001-lm.csv","./VSD001-lm.csv",method = "w")
#' download.file(url="https://github.com/marcelluethi/
#' statismo-shaperegistration/raw/master/data/VSD002-lm.csv","./VSD002-lm.csv",method = "w")
#' ref <- read.vtk("VSD001_femur.vtk")
#' tar <- read.vtk("VSD002_femur.vtk")

#' ref.lm <- as.matrix(read.csv("VSD001-lm.csv",row.names=1,header = FALSE))
#' tar.lm <- as.matrix(read.csv("VSD002-lm.csv",row.names=1,header = FALSE))
#' Kernels <- SumKernels(GaussianKernel(50,50),IsoKernel(0.1,ref))
#' mymod <- statismoModelFromRepresenter(ref,kernel=Kernels)
#' mymodC <- RvtkStatismo::statismoConstrainModel(mymod,tar.lm,ref.lm,2)
#' fit <- modelFitting(mymodC,tar,iterations = 15)
#' #or without landmarks but instead with some icp steps
#' taricp <- icp(tar,ref,iterations = 50,type="s",getTransform = TRUE)
#' taricpAff <- icp(taricp$mesh,ref,iterations = 50,type="a",getTransform = TRUE)
#' ##get affine transform
#' combotrafo <- taricpAff$transform%*%taricp$transform
#' fit2 <- modelFitting(mymod,taricpAff$mesh,iterations = 15)
#' ## revert affine transforms
#' fit2aff <- applyTransform(fit2$mesh,combotrafo,inverse=TRUE)
#' }
#' @details this function fits a statistical model to a target mesh using a l-bfgs optimizer to minimize the symmetric mean squared distance between meshes.
#' @note needs RvtkStatismo installed
#' @importFrom lbfgs lbfgs
#' @export
modelFitting <- function(model, tarmesh, iterations=5,lbfgs.iter=5,lbfgs.orthantwise_c=0, symmetric=c(0,1,2),refdist=1e5,tardist=1e5,rho=pi/2,sdmax=NULL,mahaprob=c("none","chisq","dist"),initparams=NULL,k=50,threads=0,method="lbfgs",silent=FALSE,...) {
    if (!requireNamespace("RvtkStatismo"))
        stop("you need to install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
    if (!is.null(initparams)) {
        if (length(initparams) == length(RvtkStatismo::GetPCAVarianceVector(model)))
            vars <- initparams
        else {
            warning(paste("length of initparams and number model parameters differ"))
            vars <- rep(0,length(RvtkStatismo::GetPCAVarianceVector(model)))
        }
    } else {
        vars <- rep(0,length(RvtkStatismo::GetPCAVarianceVector(model)))
    }
    if (length(lbfgs.orthantwise_c) == 1)
        lbfgs.orthantwise_c <- rep(lbfgs.orthantwise_c,iterations)
    else if (length(lbfgs.orthantwise_c) != iterations)
        stop("lbfgs.orthantwise_c must be of the same length as iterations")
    
    Aorig <- RvtkStatismo::GetPCABasisMatrix(model)
    A <- clost <- NULL
    B <- tarclost <- NULL
    mv <- RvtkStatismo::GetMeanVector(model)
    refind <- ((1:(length(mv)/3)) -1 )*3
    refind <- cbind(refind+1,refind+2,refind+3)
    symmetric <- symmetric[1]
    tarmesh <- vcgUpdateNormals(tarmesh)
    if (! symmetric %in% c(0:2))
        stop("symmetric must be 0,1 or 2")
    for ( i in 1:iterations) {
        ## to target
        mm <- vcgUpdateNormals(RvtkStatismo::DrawSample(model,vars))
        if (symmetric %in% c(0,1)) {
            cc <- vcgClostKD(mm,tarmesh,sign = FALSE,angdev=rho,k=k,threads=threads)
            ncref <- as.logical(normcheck(cc,mm) < rho)
            distgoodref <- as.logical(abs(cc$quality) <= refdist)
            goodref <- sort(which(as.logical(distgoodref*ncref)))
            refindtmp <- as.vector(t(refind[goodref,]))
            clost <- as.vector(cc$vb[1:3,goodref])-mv[refindtmp]
            A <- Aorig[refindtmp,]
        }
            ## from target
        if (symmetric %in% c(0,2)) {
            tarGet <- vcgKDtree(mm,tarmesh,k=1,threads=threads)
            dummynorms <- list(normals=mm$normals[,tarGet$index])
            nctar <- as.logical(normcheck(dummynorms,tarmesh) < rho)
            distgoodtar <- as.logical(abs(tarGet$distance) <= tardist)
            goodtar <- sort(which(as.logical(distgoodtar*nctar)))
            good <- sort(which(as.logical(distgoodtar*nctar)))
            ##good <- which(tarGet$distance < tardist)
            tarclostind <- tarGet$index[good]
            inds <- (tarclostind-1)*3
            inds <- cbind(inds+1,inds+2,inds+3)
            inds <- as.vector(t(inds))
            B <- Aorig[inds,]
            tarclost <- as.vector(tarmesh$vb[1:3,good])-mv[inds]
        }
        ##run optimization
        method <- match.arg(method[1],c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                      "Brent","lbfgs"))
        if (method != "lbfgs")
            out <- optim(par=vars,objectiveMSQ,objectiveMSQ.grad,A=A,clost=clost,B=B,tarclost=tarclost,method = method,...)
        else
            out <- lbfgs(objectiveMSQ,objectiveMSQ.grad,vars=vars,A=A,clost=clost,B=B,tarclost=tarclost,max_iterations = lbfgs.iter,invisible=1,orthantwise_c=lbfgs.orthantwise_c[i],...)
        vars <- out$par
        if (!is.null(sdmax)) {
            vars <-  constrainParams(vars,sdmax=sdmax,mahaprob = mahaprob)
        }
        
        
        if (!silent) {
            cat(paste("iteration", i,"\n"))
            cat(paste("MSE between correspondences:",objectiveMSQ(vars,clost,A,B,tarclost),"\n"))
        }
    }
    estim <- RvtkStatismo::DrawSample(model,vars)
    return(list(mesh=estim,par=vars))
}

constrainParams <- function(alpha,sdmax=3,mahaprob=c("none","chisq","dist")) {
    
    mahaprob <- match.arg(mahaprob,c("none","chisq","dist"))
    if (mahaprob != "none") {
        sdl <- length(alpha)
        probs <- sum(alpha^2)
        if (mahaprob == "chisq") {
            Mt <- qchisq(1-2*pnorm(sdmax,lower.tail=F),df=sdl)
            probs <- sum(alpha^2)
        } else if (mahaprob == "dist") {
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
    return(alpha)
}


#' minimize mean squared distance between model and a point-cloud with correspondences
#'
#' minimize mean squared distance between model and a point-cloud with correspondences
#' @param clost matrix or mesh3d
#' @param model statismo model of class pPCA
#' @param iterations integer: max number of iterations passed to lbfgs
#' @param initpar initial estimate of the model parameters
#' @param use integer vector: which points to use
#' @param sdmax constrain parameters (normalized PC-scores) to be within +- sdmax
#' @param mahaprob character: if != "none", use mahalanobis-distance to determine overall probability (of the shape projected into the model space."chisq" uses the Chi-Square distribution of the squared Mahalanobisdistance, while "dist" restricts the values to be within a multi-dimensional sphere of radius \code{sdmax}. If FALSE the probability will be determined per PC separately.
#' @param ... additional parameters to be passed to \code{\link{lbfgs}}.
#' @return
#' \item{par}{the model's parameters}
#' \item{mesh}{the fitted mesh}
#' @export
miniSQmodel <- function(clost,model,iterations=10,initpar=NULL,use=NULL,sdmax=NULL,mahaprob=c("none","chisq","dist"),...) {
     if (!requireNamespace("RvtkStatismo"))
         stop("you need to install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
    if (is.null(initpar))   
        vars <- rep(0,length(RvtkStatismo::GetPCAVarianceVector(model)))
    else
        vars <- initpar
     Aorig <- RvtkStatismo::GetPCABasisMatrix(model)
     mv <- RvtkStatismo::GetMeanVector(model)
     if (inherits(clost,"mesh3d"))
         clost <- vert2points(clost)
     A <- B <- tarclost <- NULL
     if (is.null(use))
         use <- 1:nrow(clost)
     refind <- ((1:(length(mv)/3)) -1 )*3
     refind <- cbind(refind+1,refind+2,refind+3)
     refindtmp <- as.vector(t(refind[use,]))
     clost <- as.vector(t(clost[use,]))-mv[refindtmp]
     A <- Aorig[refindtmp,]
     out <- lbfgs(objectiveMSQ,objectiveMSQ.grad,vars=vars,A=A,clost=clost,B=B,tarclost=tarclost,max_iterations = iterations,invisible=1,...)
     vars <- out$par
     if (!is.null(sdmax))
         vars <-  constrainParams(vars,sdmax=sdmax,mahaprob = mahaprob)
         
     estim <- RvtkStatismo::DrawSample(model,vars)
     return(list(mesh=estim,par=vars))
 }



