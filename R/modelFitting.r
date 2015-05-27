objectiveMSQ <- function(x,clost,A,B,tarclost) {
    Ax <- A%*%x
    out <- sum((Ax-clost)^2)/(nrow(Ax)/3)
    if (!is.null(B)) {
        Bx <- B%*%x
        out <- out+sum((Bx-tarclost)^2)/(nrow(Bx)/3)
    }
    return(out)
}


objectiveMSQ.grad <- function(x,clost,A,B,tarclost) {
    Ax <- A%*%x
    Axb <- Ax-clost
    grad <- (2*t(A)%*%Axb)/(nrow(Ax)/3)
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
#' @param tardist maximal distance from target to model to be considered
#' @param symmetric logical: if TRUE, the symmetric distance will be minimized.
#' @param silent logical: if TRUE output will be suppressed.
#' @return
#' \item{par}{the model's parameters}
#' \item{mesh}{the fitted mesh}
#' @examples
#' require(RvtkStatismo)
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD001_femur.vtk","./VSD001_femur.vtk",method = "w")
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD002_femur.vtk","./VSD002_femur.vtk",method = "w")
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD001-lm.csv","./VSD001-lm.csv",method = "w")
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD002-lm.csv","./VSD002-lm.csv",method = "w")
#' 
#' ref <- read.vtk("VSD001_femur.vtk")
#' tar <- read.vtk("VSD002_femur.vtk")
#' ref.lm <- as.matrix(read.csv("VSD001-lm.csv",row.names=1))
#' tar.lm <- as.matrix(read.csv("VSD002-lm.csv",row.names=1))
#' mymod <- statismoModelFromRepresenter(ref,kernel=list(c(50,50)),ncomp = 100,isoScale = 0.1)
#' mymodC <- statismoConstrainModel(mymod,tar.lm,ref.lm,2)
#' fit <- modelFitting(mymodC,tar,iterations = 15)
#' @details this function fits a statistical model to a target mesh using a l-bfgs optimizer to minimize the symmetric mean squared distance between meshes.
#' @note needs RvtkStatismo installed
#' @importFrom lbfgs lbfgs
#' @export
modelFitting <- function(model, tarmesh, iterations=5,symmetric=TRUE,tardist=1e5,silent=FALSE) {
    if (!require(RvtkStatismo))
        stop("you need to install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
    vars <- rep(0,length(GetPCAVarianceVector(model)))
    A <- GetPCABasisMatrix(model)
    B <- tarclost <- NULL
    mv <- GetMeanVector(model)
    for ( i in 1:iterations) {
        ## to target
        mm <- DrawSample(model,vars)
        cc <- vcgClostKD(mm,tarmesh,sign = FALSE)
        clost <- as.vector(cc$vb[1:3,])-mv
        ## from target
        if (symmetric) {
            tarGet <- vcgKDtree(mm,tarmesh,k=1)
            good <- which(tarGet$distance < tardist)
            tarclostind <- tarGet$index[good]
            inds <- (tarclostind-1)*3
            inds <- cbind(inds+1,inds+2,inds+3)
            inds <- as.vector(t(inds))
            B <- A[inds,]
            tarclost <- as.vector(tarmesh$vb[1:3,good])-mv[inds]
        }
        ##run optimization
        out <- lbfgs(objectiveMSQ,objectiveMSQ.grad,vars=vars,A=A,clost=clost,B=B,tarclost=tarclost,max_iterations = 5,invisible=1)
        vars <- out$par
        if (!silent) {
            cat(paste("iteration", i,"\n"))
            cat(paste("MSQ:",objectiveMSQ(vars,clost,A,B,tarclost),"\n"))
        }
    }
    estim <- DrawSample(model,vars)
    return(list(mesh=estim,par=vars))
}
