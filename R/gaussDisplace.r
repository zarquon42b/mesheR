## @export gaussDisplace
#' @importFrom Rvcg vcgUpdateNormals
gaussDisplace <- function(mesh1,mesh2,sigma,gamma=2,W0,f,oneway=F,nh=NULL,tol=0,pro=c("morpho","vcg","kd"),k0=50,prometh=1,rhotol=NULL,border=FALSE,horiz.disp=NULL,angclost=FALSE,bbox=NULL,threads=1,smoothtype=0,...) {
### the workhorse function running in each iteration of gaussDisplMesh3d
    ## set projection function according to input request
    pro <- substring(pro[1],1L,1L)
    if (pro == "k") {
        project3d <- vcgClostKD
    } else if (pro == "m") {
        protmp <- function(x,y,sign=F,...) {
            out <- closemeshKD(x,y,sign=sign,method=prometh,...)
            return(out)
        }
        project3d <- protmp
    } else if (pro =="v")
        project3d <- vcgClost
    
    angdev <- ifelse((is.null(rhotol) || !angclost),0,rhotol)
    
    rc <- 0
    out <- NULL
    t0 <- Sys.time()
    
    M0 <- t(mesh2$vb[1:3,])
    S0 <- t(mesh1$vb[1:3,])
    Spro <- project3d(mesh1,mesh2,sign=F,angdev=angdev,k=k0,tol=tol,borderchk=!border,threads=threads)
    
    S <- vert2points(Spro)
    ## get symmetric distances and displacement field between meshes
    if (oneway) {
        M <- vert2points(mesh2)
    } else {
        Mpro <- project3d(mesh2,mesh1,sign=F,angdev=angdev,k=k0,tol=tol,borderchk=!border,threads=threads)
        M <- vert2points(Mpro)
    }
    ## get neighbourhood for each point to minimize calculation time
    if (!is.null (nh)) {
        clostW <- vcgKDtree(S,W0,k=nh,threads=threads)
        clostIndW <- clostW$index-1L
        clostIndWdistance <- clostW$distance
        if (!oneway) {
            clostP <- vcgKDtree(M,W0,k=nh,threads=threads)
            clostIndP <- clostP$index-1L
            clostIndPdistance <- clostP$distance
        } else {
            clostIndPdistance <- clostIndP <- matrix(0,dim(W0)[1],nh)
        }
    }
    rt0 <- rep(0,nrow(S))
    rt1 <- rep(0,nrow(M))
    if (!is.null(rhotol)) {
        rc <- rhotol
        rt0 <- normcheck(mesh1,Spro)
        if (!oneway)
            rt1 <- normcheck(mesh2,Mpro)
    }
    if (!is.null(horiz.disp)) {
        if (is.null(rhotol))
            rc <- horiz.disp
        tmp <- list();tmp$normals <- mesh1$vb[1:3,]-Spro$vb[1:3,]
        hordev0 <- normcheck(mesh1,tmp)
        rt0[which(hordev0 > horiz.disp)] <- 4
        if (!oneway) {
            tmp <- list();tmp$normals <- mesh2$vb[1:3,]-Mpro$vb[1:3,]
            hordev1 <- normcheck(mesh2,tmp)
            rt1[which(hordev1 > horiz.disp)] <- 4
        }
    }
    if (!is.null(bbox)) {
        badrange <- outsideBBox(mesh1,bbox)
        
        rt0[badrange] <- 4
        
        if (!oneway) {
            badrange <- outsideBBox(Mpro,bbox)
            rt1[badrange] <- 4
        }
    }
    t3 <- Sys.time()
    D1 <- S-S0
    D2 <- M-M0
    if (!border) {
        if (is.null(rhotol))
            rc <- pi
        if (pro %in% c("v","k")) {
            rt0[as.logical(Spro$border)] <- 4
            if (!oneway)
                rt1[as.logical(Mpro$border)] <- 4
        } else {
            bordtmp <- vcgBorder(mesh2)
            rt0[which(Spro$faceptr %in% which(as.logical(bordtmp$borderit)))] <- 4
            if (!oneway) {
                bordtmp <- vcgBorder(mesh1)
                rt1[which(Mpro$faceptr %in% which(as.logical(bordtmp$borderit)))] <- 4
            }
        }
    }
    tol <- tol^2
### make multicore 
    if (nh < 1)
        stop ("neighbourhood must be at least 1")
    out <- .Call("displaceGauss",W0,M0,D1,D2,sigma,gamma,clostIndW,clostIndWdistance,clostIndP,clostIndPdistance,tol=tol,rt0,rt1,rc,oneway,smoothtype,threads,PACKAGE="mesheR")
    addit <- W0+out
    return(list(addit=addit))
}



#' map two surface meshes using smoothed displacement fields
#' 
#' Map a reference mesh onto a target surface using displacement fields.
#' 
#' This function implements the mesh matching method suggested by Moshfeghi et
#' al. and Bryan et al.. Additional mechanisms for controlling and restricting
#' the displacement smoothing are implemented
#' 
#' @param x reference mesh: triangular mesh of class "mesh3d"or of class BayesDeform created by createBayes to restrict based on a known distribution. To use this option the package RvtkStatismo \url{https://github.com/zarquon42b/RvtkStatismo} has to be installed. If x is a model, it works best if mesh2 is already aligned to the model's mean.
#' @param mesh2 An object of class mesh3d used as target mesh. Mesh resolution
#' should be ~1.5.
#' @param iterations Iterations of displacement. Default is 10.
#' @param smooth Integer: smooth the resulting mesh at the \code{smooth-th} iteration. E.g. a value of 2 will smooth the resulting mesh at each 2nd iteration. Not to be confused with \code{displacementsmooth}!! Default is NULL, no smoothing.
#' @param smoothit integer: smoothing steps.
#' @param smoothtype Type of smoothing: Taubin, Laplacian, or HClaplacian. For
#' details see \code{vcgSmooth}.
#' @param sigma starting value of kernel bandwidth/B-spline support. For all kernels except B-spline, sigma controls the importance of the neighbourhood by defining the bandwidth of the smoothing kernel. For B-spline it defines the support (the higher, the "wobblier" the deformation field can become.
#' @param displacementsmooth kernel function for smoothing are "Gauss","Laplace", "Exponential" and "Bspline" (or any abbreviation thereof).
#' @param gamma dampening factor controlling displacement strength. The smoothed displacement vector for each vertex is divided by \code{gamma}. The larger \code{gamma}, the slower the approximation.
#' @param f parameter controlling iterative decrease (if > 1, increase else) of \code{sigma} making the displacement locally more elastic with each iteration.
#' Starting with \code{sigma}, this parameter for the k-th iteration is \code{sigma *f ^(-k)}
#' @param oneway logical: only displace towards the target without taking into
#' account the displacement from the target.
#' 
#' @param lm1 A k x 3 matrix containing landmarks corrresponding to mesh1 for
#' initial rotation of mesh1 onto mesh2.
#' @param lm2 A k x 3 matrix containing landmarks corrresponding to mesh2 for
#' initial rotation of mesh1 onto mesh2.
#' @param rigid named list. Passing parameters to \code{\link{icp}}, for rigid registration. If landmarks are provided and only those should count, set rigid$iterations=0.
#' @param similarity named list. Passing parameters to \code{\link{icp}}, for similarity registration (rigid +scaling). If landmarks are provided and only those should count, set similarity$iterations=0 (and rigid=NULL).
#'@param affine named list. Passing parameters to \code{\link{icp}}, for affine registration. If landmarks are provided and only those should count, set similarity$iterations=0 (with rigid=NULL and similarity=NULL)
#' @param tps logical: if TRUE and landmarks are provided, the reference will be mapped to the target using a Thin-Plate Spline interpolation. Overrides \code{rigid},\code{affine} and \code{similarity}.
#' @param nh Integer: neighbourhood (number vertices) for controlling
#' displacement smoothing, default is 150/mesh resolution.
#' @param toldist Integer: Exclude everything from the whole procedure with a
#' greater distance from initial point than toldist. 0 disables this feature.
#' @param pro which projection method to use: "kd" = \code{\link{vcgClostKD}}, "m"= \code{\link{closemeshKD}}
#' from Morpho; "v"= \code{\link{vcgClost}} from package Rvcg.
#' @param k0 Integer: argument passed to closemeshKD (will be argument "k" in
#' \code{\link{closemeshKD}} .
#' @param prometh argument passed to closemeshKD.  Integer: 0 or 1. If
#' prometh=0, take closest point for displacement. If prometh=1, do not just
#' take the closest point, but for two absolut distances which are the same,
#' take the point which is orthogonal to the closest face see Moshfeghi 1994).
#' @param angtol numeric: If the angle between hit points' normals and the
#' starting points' normals exceeds this threshold the displacement vector will
#' be discarded. 
#' @param border Logical: if TRUE, displacement vectors hitting mesh borders
#' are discarded.
#' @param horiz.disp numeric: If the angle between hit points' normals
#' (independent of its orientation) and the distance vector between hit point
#' and starting points exceeds this threshold, the displacement vector will be
#' discarded. Reduces distortion especially at mesh borders.
#' @param useiter logical: if AmbergK and AmbergLambda are set and useiter is TRUE, the minimization will be performed using the latest iteration (slower). The orginal reference otherwise.
#' @param AmbergK a single integer or an integer vector vector containing the \code{k0}-value (normal slackness) for each iteration for a smooth Deformation using \code{\link{AmbergDeformSpam}}.
#' @param AmbergLambda as single numeric value or a numeric vector containing the \code{lambda}-value for each iteration for a smooth Deformation using \code{\link{AmbergDeformSpam}}.
#' @param tol convergence threshold: if RMSE between iterations is below tol, the function stops.
#' @param useConstrained logical: if TRUE and Bayes and landmarks are defined, the landmarks are not only used to get a suitable reference but the model will also be constrained by the landmarks to subsequently restrict the shape variability. If FALSE, the full model is used.
#' @param angclost if TRUE, the closest k faces will be evaluated and the closest with the appropriate normal angle will be selected.
#' @param noinc after each iteration the RMSE between target and moving image is calculated and if this value increases compared to a previous value, the matching stops. Can be useful when matching a statistical model to a partial shape.
#' @param silent logical suppress messages
#' @param visualize logical: if TRUE the matching process is visualized
#' @param folder character: if specified, each a screenshot of each deformation step will be saved as a png file in this folder.
#' @param alpha numeric between 0 and 1 controls opacity of target mesh if visualize=TRUE.
#' @param col1 color of fix mesh (if visualize = TRUE)
#' @param col2 color of moving mesh (if visualize = TRUE)
#' @param bbox extend of the margins around the target shape to be considered.
#' @param bboxCrop extend of the bounding box around mesh1 (after alignmend) that will be cropped from target to speed things up.
#' @param threads integer: threads to use in multithreaded routines.
#' @param cb optional: callback function that takes arguments i="current iteration", distance= "distance from target to current estimate" and t.dist="average vertex displacement to last iteration"
#' @param \dots Further arguments passed to \code{nn2}.
#'
#' @return If a patch is specified:
#'  \item{mesh}{matched mesh}
#'  \item{patch}{displaced patch as specified in input.}
#' else a mesh of class "mesh3d" is returned.
#' @note Based on the closest points (constrained by the various options), the displacement field will be smoothed using kernel functions of the k-closest displacement vectors.
#'  The smoothing kernels are  "Gauss","Laplace", "Exponential" and "Bspline". The displacement at point \code{x} will be the weighted displacment vectors of the k-closest displacement vectors. Be \code{d} the distance to a neightbouring point, the weight will be calculated as:
#' 
#' Gaussian:  \eqn{w(d) = exp(\frac{-d^2}{2\sigma^2})}{w(d) = exp(-d^2/2*sigma^2)}
#' 
#' Laplacian: \eqn{w(d) = exp(\frac{-d}{\sigma})}{w(d) = exp(-d/sigma)}
#'
#' Exponential: \eqn{w(d) = exp(\frac{-d}{2\sigma^2})}{w(d) = exp(-d/2*sigma^2)}
#' @author Stefan Schlager
#' @seealso \code{\link{meshres}}, \code{\link{vcgClost}},
#' \code{\link{vcgBorder}}, \code{\link{icp}}, \code{\link{vcgSmooth}}
#' @references Bryan, R., Mohan, P. S., Hopkins, A., Galloway, F., Taylor, M.,
#' and Nair, P. B. 2010. Statistical modelling of the whole human femur
#' incorporating geometric and material properties. Medical Engineering &amp;
#' Physics, 32(1):57-65.
#' 
#' Moshfeghi, M., Ranganath, S., and Nawyn, K. 1994. Three-dimensional elastic
#' matching of volumes. IEEE Transactions on Image Processing: A Publication of
#' the IEEE Signal Processing Society, 3(2):128-138.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' require(Morpho)
#' data(nose)##load data
#' ##warp a mesh onto another landmark configuration:
#' longnose.mesh <- tps3d(shortnose.mesh,shortnose.lm,longnose.lm)
#' ### result won't be too good as the surfaces do stronly differ.
#' ## we start with an affine transformation initiated by landmarks
#' affine <- list(iterations=20,subsample=100,rhotol=pi/2,uprange=0.9)
#' match <- gaussMatch(shortnose.mesh,longnose.mesh,lm1=shortnose.lm,
#'                    lm2=longnose.lm,gamma=2,iterations=10,smooth=1,smoothtype="h",
#'                    smoothit=10,nh=50,angtol=pi/2,affine=affine,sigma=20)
#'
#' \dontrun{
#' ## now a more cautiously approach using a larger neighbourhood
#' ## (400 instead of 50) and no intermediary smoothing:
#' matchNoSmooth <- gaussMatch(shortnose.mesh,longnose.mesh,lm1=shortnose.lm,
#'                    lm2=longnose.lm,gamma=2,iterations=20,,nh=400,angtol=pi/2,
#'                    affine=affine,sigma=20)
#' }
#' @importFrom Rvcg vcgClostKD vcgKDtree vcgMeshres
#' @importFrom rgl rgl.ids
#' @seealso \code{\link{outsideBBox},\link{getMeshBox} }
#' @export
#'
#' @useDynLib mesheR
gaussMatch <- function(x,mesh2,iterations=10,smooth=NULL,smoothit=10,smoothtype=c("taubin","laplace","HClaplace"),sigma=20,displacementsmooth=c("Gauss","Laplace","Exponential"),gamma=2,f=1.2,oneway=F,lm1=NULL,lm2=NULL,rigid=NULL, similarity=NULL, affine=NULL,tps=FALSE,nh=NULL,toldist=0,pro=c("kd","vcg","morpho"),k0=50,prometh=1,angtol=pi/2,border=FALSE,horiz.disp=NULL,useiter=FALSE,AmbergK=NULL,AmbergLambda=NULL,tol=1e-5, useConstrained=TRUE, angclost=TRUE,noinc=FALSE,silent=FALSE, visualize=FALSE,folder=NULL,alpha=0.7,col1="red",col2="white",bbox=NULL,bboxCrop=NULL,threads=0,cb=NULL,...) {
    typeargs <- c("gauss","laplace","exponential","bspline")
    displacementsmooth <- match.arg(tolower(displacementsmooth[1]),typeargs)
    displacementsmooth <- match(displacementsmooth,typeargs)-1
    sigma0 <- sigma
    if (inherits(x, "mesh3d")) {
        mesh1 <- x
        Bayes <- NULL
    } else if (inherits(x, "BayesDeform"))
        Bayes <- x
    else
        stop("x must be an object of class mesh3d or BayesDeform")
    if (!is.null(Bayes)) {
        if (!requireNamespace("RvtkStatismo"))
            stop("for using the option Bayes, please install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
        mesh1 <- RvtkStatismo::DrawMean(Bayes$model)
    }
    
    
    if (!is.null(angtol)) {
        mesh1 <- vcgUpdateNormals(mesh1)
        mesh2 <- vcgUpdateNormals(mesh2)
    }
    Amberg <- ambergsingle <- FALSE
    ##setup variables
    if (!is.null(AmbergK) && !is.null(AmbergLambda)) {
        AmbergK <- round(AmbergK)# make sure k is integer - otherwise RAM overkill
        ambergsingle <- FALSE
        if (length(AmbergK) == 1) {
            AmbergK <- rep(AmbergK,iterations)
            ambergsingle <- TRUE
        } else if (length(AmbergK) != iterations)
            stop("AmbergK must be vector of length 'iterations'")
        
        if (length(AmbergLambda) == 1)
            AmbergLambda <- rep(AmbergLambda,iterations)
        else if (length(AmbergLambda) != iterations)
            stop("AmbergLambda must be vector of length 'iterations'")
        else
            ambergsingle <- FALSE
        Amberg <- TRUE
        Hchol <- NULL
    }
    ## clean input mesh
    if(length(unrefVertex(mesh1)) > 0 )
        mesh1 <- rmUnrefVertex(mesh1)
    
    if (is.null(nh)) {
        nh=ceiling(150/vcgMeshres(mesh1)$res)
        if (!silent)
            cat(paste("\nneighbourhood is set to",nh,"\n***************\n"))
    }
    t.dist <- 1e12
    hasLM <- FALSE
    if (!is.null(lm1) && !is.null(lm2)) {
        hasLM <- TRUE
        bary <- vcgClost(lm1,mesh1,barycentric = T)
    }
    if (!is.null(Bayes$initparams)) {
        mesh1 <- vcgUpdateNormals(RvtkStatismo::DrawSample(Bayes$model,Bayes$initparams))
        if (hasLM)
            lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
                                        #wire3d(mesh1);spheres3d(lm1)
                                        #return(1)
    }
    
    ## do icp matching
    lmModel <- NULL
    if (hasLM) {   ## case: landmarks are provided
        
        if (!is.null(Bayes)) {
            if (Bayes$align)
                lm2tmp <- rotonto(lm1,lm2,scale=Bayes$model@scale,reflection=FALSE)$yrot
            else
                lm2tmp <- lm2
            if (!silent)
                cat("constraining model\n")
            constMod <- RvtkStatismo::statismoConstrainModel(Bayes$model,lm2tmp,lm1,Bayes$ptValueNoise)
            if (useConstrained) {
                Bayes$model <- constMod
                if (is.null(Bayes$initparams)) {
                    mesh1 <- vcgUpdateNormals(RvtkStatismo::DrawMean(Bayes$model))
                    lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
                }
            }
        }
        if (tps) {
            mesh1 <- tps3d(mesh1,lm1,lm2,threads=threads)
            lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
        } else {
            if (is.null(rigid) && is.null(affine) && is.null(similarity)) {
                if (is.null(Bayes))                
                    rigid <- list(iterations=0)
                else if (Bayes$align)
                    rigid <- list(iterations=0)
            }
            if (!is.null(rigid)) { ##perform rigid icp-matching
                rigid$lm1 <- lm1
                rigid$lm2 <- lm2
                mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent,threads=threads)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
            if (!is.null(similarity)) {##similarity matching
                if (is.null(rigid)) {
                    similarity$lm1 <- lm1
                    similarity$lm2 <- lm2
                }
                mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent,threads=threads)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
            if (!is.null(affine)) {##similarity matching
                if (is.null(rigid) && is.null(similarity)) {
                    affine$lm1 <- lm1
                    affine$lm2 <- lm2
                }
                mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent,threads=threads)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
        }
        if (!is.null(Bayes)) 
            mesh1 <- vcgUpdateNormals(RvtkStatismo::PredictSample(Bayes$model,mesh1,representer = T,align=Bayes$align,sdmax=Bayes$sdmax[1],mahaprob=Bayes$mahaprob))
        
                                        #mesh1 <- vcgUpdateNormals(RvtkStatismo::PredictSample(Bayes$model,mesh1,representer = T,lmDataset=lm1,lmModel=lmModel,align=TRUE))
        
    } else {
        if (!is.null(rigid) || !is.null(affine) || !is.null(similarity)) {
            if (!is.null(rigid)) ##perform rigid icp-matching
                mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent,threads = threads)
            if (!is.null(similarity))##similarity matching
                mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent,threads = threads)
            if (!is.null(affine))##similarity matching
                mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent,threads = threads)
        }
    }
    if (!is.null(bboxCrop)) {
        mesh2 <- cropOutsideBBox(mesh1,mesh2,extend=bboxCrop)
        if (!silent)
            cat("cropping target mesh\n")
    } 
    
    if (!is.null(bbox))
        bbox <- getMeshBox(mesh2,extend=bbox)
    
    if (visualize) {
        rglid <- NULL
        if (!length(rgl.ids()$id)) 
            open3d()
        else {
            rgl.bringtotop()
            rgl.clear()
        }
        bb <- meshcube(mesh1)
        bmean <- apply(bb,2,mean)
        bb <- t(((t(bb)-bmean)*2)+bmean)
        points3d(bb,col="white",alpha=0)
        shade3d(mesh2,col=col1,specular=1,alpha=alpha)
        if (!is.null(rglid))
            rgl.pop(id=rglid)
        rglid <- shade3d(mesh1,col=col2,front="lines", back="lines")
        
        if (!is.null(folder)) {
            if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") 
                folder <- paste(folder,"/",sep="")
            dir.create(folder,showWarnings=F)
            movie <- paste(folder,"deformation",sep="")
            
            npics <- nchar(iterations+1)
            ndec <- paste0("%s%0",npics,"d.png")
        }
        if (interactive())
            readline("please select viewpoint\n")
        
        
        if (!is.null(folder)) {
            filename <- sprintf("%s%04d.png", movie, 1)
            rgl.snapshot(filename,fmt="png")
            movcount <- 2
        }
    }
    if (Amberg) {
        if (!useiter)
            S <- createS(mesh1)
        else
            S <- NULL
        meshorig <- mesh1
    }
    ## elastic matching starts
    if (!silent)
        cat("starting elastic matching\n****************\n")
    i <- 1
    distance <- 1e10
    while (i <= iterations && t.dist > tol ) {
        time0 <- Sys.time()
        if (!is.null(smooth) && i > 1) {
            if (i %% smooth == 0) {
                if (!silent)
                    cat("smoothing step\n")
                mesh1 <- vcgSmooth(mesh1,type=smoothtype,iteration=smoothit)
                                        #if (!silent)
                                        #cat("smoothing finished\n")
            }
        }
        vb0 <- vert2points(mesh1)
        sigma <- sigma0*f^(-(i-1))
        if (!silent)
            cat(paste("sigma =",round(sigma,3),"\n"))
        ## call the workhorse doing the displacement
        tmp <- gaussDisplace(mesh1,mesh2,sigma=sigma,gamma=gamma,f=f,W0=vert2points(mesh1),nh=nh,tol=toldist,pro=pro,k0=k0,prometh=prometh,rhotol=angtol,border=border,oneway=oneway,horiz.disp = horiz.disp,bbox=bbox,angclost=angclost,threads=threads,smoothtype=displacementsmooth,...)
        
        tmpold <- mesh1
        if (!is.null(Bayes) && length(Bayes$sdmax) >= i) {
            if (!is.null(Bayes$wt)) {
                mesh0 <- mesh1
                mesh0$vb[1:3,] <- t(tmp$addit)
                wt <- Bayes$wt[i]
                wts <- c(1,wt)
                wts <- wts/sum(wts)
                tmpmesh <- RvtkStatismo::PredictSample(Bayes$model,lmDataset=lm1,lmModel=lmModel,dataset=mesh0,representer=TRUE, sdmax=Bayes$sdmax[i],align=Bayes$align,mahaprob=Bayes$mahaprob)
                tmp$addit <- t(wts[1]*mesh0$vb[1:3,]+wts[2]*tmpmesh$vb[1:3,])
                
            } else {
                tmp$addit <- RvtkStatismo::PredictSample(Bayes$model,tmp$addit,FALSE, sdmax=Bayes$sdmax[i],mahaprob=Bayes$mahaprob,align=Bayes$align)
            }
            
        }
        
        if (Amberg) {
            if (!useiter) 
                mytry <- try(ambtry <- AmbergDeformSpam(meshorig,vert2points(meshorig),tmp$addit,lambda=AmbergLambda[i],k0=AmbergK[i],S=S,Hchol=Hchol),FALSE)
            else
                mytry <- try(ambtry <- AmbergDeformSpam(mesh1,vert2points(mesh1),tmp$addit,lambda=AmbergLambda[i],k0=AmbergK[i],S=S,Hchol=Hchol),FALSE)

            if (!inherits(mytry,"try-error")) {
                if (ambergsingle && !useiter) 
                    Hchol <- ambtry$Hchol
                
                mesh1 <- ambtry$mesh
            }
        } else {
            mesh1$vb[1:3,] <- t(tmp$addit)
        }
        ## check if distance increases
        distance_old <- distance
        distance <- mean(vcgClostKD(mesh2,mesh1,k0=2,sign=F,threads=threads)$quality)
        if (distance > distance_old && !is.null(Bayes) && noinc) {
            cat("\n=========================================\n")
            message(paste(" Info: Distance is increasing, matching stopped after ",i,"iterations\n"))
            i <- 1e10
            mesh1 <- tmpold
        }
        mesh1 <- vcgUpdateNormals(mesh1)
        if (visualize) {
            
            if (!is.null(rglid))
                rgl.pop(id=rglid)
            rglid <- shade3d(mesh1,col=col2,front="lines", back="lines",)
            if (!is.null(folder)) {
                filename <- sprintf("%s%04d.png", movie, movcount)
                movcount <- movcount+1
                rgl.snapshot(filename,fmt="png")
            }
        }
        
        t.dist <- mean(sqrt(rowSums((vert2points(mesh1)-vb0)^2)))
        time1 <- Sys.time()
        gc()
        if (!silent && i < 1e10) {
            cat(paste("completed iteration",i, "in", round(as.numeric(time1-time0,unit="secs"),2), "seconds\n"))
            cat(paste(" Info: Average distance to target:",distance,"\n"))
            cat(paste0("average vertex displacement to last iteration = ",t.dist,"\n"))
            cat("****************\n")
            
        }
        if (i <1e10 && !is.null(cb))
            cb(i,distance,t.dist)

        i <- i+1
        
    }
    invisible(mesh1)
}


