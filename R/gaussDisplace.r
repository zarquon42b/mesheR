## @export gaussDisplace
#' @importFrom Rvcg vcgUpdateNormals
gaussDisplace <- function(mesh1,mesh2,sigma,gamma=2,W0,f,oneway=F,k=1,nh=NULL,tol=0,pro=c("morpho","vcg"),k0=50,prometh=1,rhotol=NULL,border=FALSE,horiz.disp=NULL,...)
{
### the workhorse function running in each iteration of gaussDisplMesh3d
    ## set projection function according to input request
    pro <- substring(pro[1],1L,1L)
    if (pro == "v") {
        project3d <- vcgClostKD
    } else if (pro == "m") {
        protmp <- function(x,y,sign=F) {
            out <- closemeshKD(x,y,k=k0,sign=sign,method=prometh)
            return(out)
        }
        project3d <- protmp
    }
    rc <- 0
    out <- NULL
    t0 <- Sys.time()
    sigma0 <- sigma
    M0 <- t(mesh2$vb[1:3,])
    S0 <- t(mesh1$vb[1:3,])
    sigma <- (sigma0*f^(-k))^2
    Spro <- project3d(mesh1,mesh2,sign=F)
    S <- vert2points(Spro)
    ## get symmetric distances and displacement field between meshes
    if (oneway) {
        M <- vert2points(mesh2)
    } else {
        Mpro <- project3d(mesh2,mesh1,sign=F)
        M <- vert2points(Mpro)
    }
    ## get neighbourhood for each point to minimize calculation time
    if (!is.null (nh)) {
        clostIndW <- vcgKDtree(S,W0,k=nh)$index-1
        if (!oneway)
            clostIndP <- vcgKDtree(M,W0,k=nh)$index-1
        else
            clostIndP <- matrix(0,dim(W0)[1],nh)
    }
    rt0 <- rep(0,dim(S)[1])
    rt1 <- rep(0,dim(M)[1])
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
        hordev0 <- normcheck(mesh1,tmp,circle=FALSE)
        rt0[which(hordev0 > horiz.disp)] <- 4
        if (!oneway) {
            tmp <- list();tmp$normals <- mesh2$vb[1:3,]-Mpro$vb[1:3,]
            hordev1 <- normcheck(mesh2,tmp,circle=FALSE)
            rt1[which(hordev1 > horiz.disp)] <- 4
        }
    }
    t3 <- Sys.time()
    D1 <- S-S0
    D2 <- M-M0
    if (!border) {
        if (is.null(rhotol))
            rc <- pi
        if (pro=="v") {
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
    storage.mode(rt0) <- "double"
    storage.mode(rt1) <- "double"
    storage.mode(clostIndW) <- "integer"
    storage.mode(clostIndP) <- "integer"
    storage.mode(S0) <- "double"
    storage.mode(M) <- "double"
    storage.mode(D1) <- "double"
    storage.mode(D2) <- "double"
    storage.mode(nh) <- "integer"
    
    tol <- tol^2
### make multicore 
    
    out <- .Call("displaceGauss",W0,S0,M,D1,D2,sigma,gamma,clostIndW,clostIndP,tol=tol,rt0,rt1,rc,oneway,PACKAGE="mesheR")
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
#' @param mesh1 x reference mesh: triangular mesh of class "mesh3d"or of class BayesDeform created by createBayes to restrict based on a known distribution. To use this option the package RvtkStatismo \url{https://github.com/zarquon42b/RvtkStatismo} has to be installed. Mesh resolution should be ~1.5.
#' @param mesh2 An object of class mesh3d used as target mesh. Mesh resolution
#' should be ~1.5.
#' @param iterations Iterations of displacement. Default is 10.
#' @param smooth Integer: smoothing factor. Default is NULL, no smoothing.
#' @param smoothit integer: smoothing steps.
#' @param smoothtype Type of smoothing: Taubin, Laplacian, or HClaplacian. For
#' details see \code{\link{vcgSmooth}}
#' @param sigma starting parameter for smoothed displacement (see Moshfeghi
#' 1994). Sigma controls the importance of the neighbourhood by defining the standard-deviation for the gaussian smoothing
#' @param gamma stiffness factor controlling displacement strength. The smoothed displacement vector for each vertex is divided by \code{gamma}. The larger \code{gamma}, the slower the approximation.
#' @param f parameter controlling iterative decrease of \code{sigma} making the displacement locally more elastic with each iteration.
#' (Moshfeghi 1994). Starting with \code{sigma}, this parameter for the k-th iteration is \code{sigma *f ^(-k)}
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
#' @param nh Integer: neighbourhood (number vertices) for controlling
#' displacement smoothing, default is 150/mesh resolution.
#' @param toldist Integer: Exclude everything from the whole procedure with a
#' greater distance from initial point than toldist. 0 disables this feature.
#' @param pro which projection method to use: "m"= \code{\link{closemeshKD}}
#' from Morpho; "v"= \code{\link{vcgClost}} from package Rvcg
#' @param k0 Integer: argument passed to closemeshKD (will be argument "k" in
#' \code{\link{closemeshKD}} .
#' @param prometh argument passed to closemeshKD.  Integer: 0 or 1. If
#' prometh=0, take closest point for displacement. If prometh=1, do not just
#' take the closest point, but for two absolut distances which are the same,
#' take the point which is orthogonal to the closest face see Moshfeghi 1994).
#' @param angtol numeric: If the angle between hit points' normals and the
#' starting points' normals exceeds this threshold the displacement vector will
#' be discarded.
#' 
#' @param border Logical: if TRUE, displacement vectors hitting mesh borders
#' are discarded.
#' @param horiz.disp numeric: If the angle between hit points' normals
#' (independent of its orientation) and the distance vector between hit point
#' and starting points exceeds this threshold, the displacement vector will be
#' discarded. Reduces distortion especially at mesh borders.
#' @param AmbergK a single integer or an integer vector vector containing the \code{k0}-value (normal slackness) for each iteration for a smooth Deformation using \code{\link{AmbergDeformSpam}}.
#'  @param AmbergLambda as single numeric value or a numeric vector containing the \code{lambda}-value for each iteration for a smooth Deformation using \code{\link{AmbergDeformSpam}}.
#' @param silent logical suppress messages
#' @param Bayes optional: object of class BayesDeform created by createBayes to restrict based on a known distribution
#' @param useConstrained logical: if TRUE and Bayes and landmarks are defined, the landmarks are not only used to get a suitable reference but the model will also be constrained by the landmarks to subsequently restrict the shape variability. If FALSE, the full model is used.
#' @param \dots Further arguments passed to \code{nn2}.
#'
#' @return If a patch is specified:
#'  \item{mesh}{matched mesh}
#'  \item{patch}{displaced patch as specified in input.}
#' else a mesh of class "mesh3d" is returned.
#'
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
#' warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
#' ### result won't be too good as the surfaces do stronly differ.
#' ## we start with an affine transformation initiated by landmarks
#' affine <- list(iterations=200,subsample=100,rhotol=pi/2,uprange=0.9)
#'  match <- gaussMatch(shortnose.mesh,warpnose.long,lm1=shortnose.lm,lm2=longnose.lm,gamma=4,iterations=10,smooth=1,smoothtype="h",smoothit=10,nh=50,angtol=pi/2,affine=affine,sigma=100)
#' @importFrom Rvcg vcgClostKD vcgKDtree
#' @export
#'
#' @useDynLib mesheR
gaussMatch <- function(x,mesh2,iterations=10,smooth=NULL,smoothit=10,smoothtype=c("taubin","laplace","HClaplace"),sigma=20,gamma=2,f=1.2,oneway=F,lm1=NULL,lm2=NULL,rigid=NULL, similarity=NULL, affine=NULL,nh=NULL,toldist=0,pro=c("vcg","morpho"),k0=50,prometh=1,angtol=NULL,border=FALSE,horiz.disp=NULL,AmbergK=NULL,AmbergLambda=NULL,silent=FALSE, Bayes=NULL,useConstrained=TRUE,visualize=FALSE,folder=NULL,...)
    {
        if (inherits(x, "mesh3d")) {
            mesh1 <- x
            Bayes <- NULL
        } else if (inherits(x, "BayesDeform"))
            Bayes <- x
        else
            stop("x must be an object of class mesh3d or BayesDeform")
        if (!is.null(Bayes)) {
            if (!require(RvtkStatismo))
                stop("for using the option Bayes, please install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
            mesh1 <- DrawMean(Bayes$model)
        }
        if (!is.null(angtol)) {
            mesh1 <- vcgUpdateNormals(mesh1)
            mesh2 <- vcgUpdateNormals(mesh2)
        }
        Amberg <- FALSE
        ##setup variables
        if (!is.null(AmbergK) && !is.null(AmbergLambda)) {
            AmbergK <- round(AmbergK)# make sure k is integer - otherwise RAM overkill
            if (length(AmbergK) == 1)
                AmbergK <- rep(AmbergK,iterations)
            else if (length(AmbergK) != iterations)
                stop("AmbergK must be vector of length 'iterations'")
            
            if (length(AmbergLambda) == 1)
                AmbergLambda <- rep(AmbergLambda,iterations)
            else if (length(AmbergLambda) != iterations)
                stop("AmbergLambda must be vector of length 'iterations'")
            Amberg <- TRUE
        }
        
        
        ## clean input mesh
        if(length(unrefVertex(mesh1)) > 0 )
            mesh1 <- rmUnrefVertex(mesh1)
        
        if (is.null(nh)) {
            nh=ceiling(150/meshres(mesh1))
            if (!silent)
                cat(paste("\nneighbourhood is set to",nh,"\n***************\n"))
        }
        ## set projection function according to input request
        pro <- substring(pro[1],1L,1L)
        if (pro == "v") {
            project3d <- vcgClostKD
        } else if (pro == "m") {
            protmp <- function(x,y,sign=F) {
                out <- closemeshKD(x,y,k=k0,sign=sign)
                return(out)
            }
            project3d <- protmp
        }
        
        ## do icp matching
        if (!is.null(lm1) && !is.null(lm2)) {   ## case: landmarks are provided
            bary <- vcgClost(lm1,mesh1,barycentric = T)
            if (!is.null(Bayes)) {
                lm2tmp <- rotonto(lm1,lm2,scale=Bayes$model@scale,reflection=FALSE)$yrot
                constMod <- statismoConstrainModel(Bayes$model,lm2tmp,lm1,Bayes$ptValueNoise)
                if (useConstrained)
                    Bayes$model <- constMod
                mesh1 <- DrawMean(statismoConstrainModel(Bayes$model,lm2tmp,lm1,Bayes$ptValueNoise))
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }           
            if (is.null(rigid) && is.null(affine) && is.null(similarity))
                rigid <- list(iterations=0)
            if (!is.null(rigid)) { ##perform rigid icp-matching
                rigid$lm1 <- lm1
                rigid$lm2 <- lm2
                mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
            if (!is.null(similarity)) {##similarity matching
                if (is.null(rigid)) {
                    similarity$lm1 <- lm1
                    similarity$lm2 <- lm2
                }
                mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
            if (!is.null(affine)) {##similarity matching
                if (is.null(rigid) && is.null(similarity)) {
                    affine$lm1 <- lm1
                    affine$lm2 <- lm2
                }
                mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
            
        } else {
            if (!is.null(rigid) || !is.null(affine) || !is.null(similarity)) {
                if (!is.null(rigid)) ##perform rigid icp-matching
                    mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent)
                if (!is.null(similarity))##similarity matching
                    mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent)
                if (!is.null(affine))##similarity matching
                    mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent)
            }
        }
        if (visualize) {
             rglid <- NULL
            open3d()
            points3d(meshcube(mesh1),col="white",alpha=0)
            shade3d(mesh2,col=2,specular=1)
            if (!is.null(rglid))
                rgl.pop(id=rglid)
            rglid <- wire3d(mesh1,col="white")
           
            if (!is.null(folder)) {
                if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") 
                    folder <- paste(folder,"/",sep="")
                dir.create(folder,showWarnings=F)
                movie <- paste(folder,"deformation",sep="")
                
                npics <- nchar(iterations+1)
                ndec <- paste0("%s%0",npics,"d.png")
            }
            readline("please select viewpoint\n")
            
            
            if (!is.null(folder)) {
                filename <- sprintf("%s%04d.png", movie, 1)
                rgl.snapshot(filename,fmt="png")
                movcount <- 2
            }
         }
        ## elastic matching starts
        if (!silent)
            cat("starting elastic matching\n****************\n")
        for (i in 1:iterations) {
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
            ## call the workhorse doing the displacement
            tmp <- gaussDisplace(mesh1,mesh2,sigma=sigma,gamma=gamma,f=f,W0=vert2points(mesh1),nh=nh,k=i,tol=toldist,pro=pro,k0=k0,prometh=prometh,rhotol=angtol,border=border,oneway=oneway,horiz.disp = horiz.disp,...)
            
            if (Amberg) {#smooth deformation
                tmpmesh <- mesh1
                tmpmesh$vb[1:3,] <- t(tmp$addit)
                tmpmesh <- vcgUpdateNormals(mesh1)
                mytry <- try(mesh1 <- AmbergDeformSpam(mesh1,vert2points(mesh1),tmp$addit,lambda=AmbergLambda[i],k0=AmbergK[i])$mesh,TRUE)
               
            } else
                mesh1$vb[1:3,] <- t(tmp$addit)


            if (!is.null(Bayes) && length(Bayes$sdmax) >= i) {
                if (!is.null(Bayes$wt)) {
                    wt <- Bayes$wt
                    wts <- c(1,wt)
                    wts <- wts/sum(wts)
                    tmpmesh <- PredictSample(Bayes$model,mesh1,TRUE, sdmax=Bayes$sdmax[i],align=TRUE)
                    mesh1$vb[1:3,] <- wts[1]*mesh1$vb[1:3,]+wts[2]*tmpmesh$vb[1:3,]
                    
                } else
                    mesh1 <- PredictSample(Bayes$model,mesh1,TRUE, sdmax=Bayes$sdmax[i],align=Bayes$align)

            }
            mesh1 <- vcgUpdateNormals(mesh1)
            if (visualize) {
                    
                    if (!is.null(rglid))
                        rgl.pop(id=rglid)
                    rglid <- wire3d(mesh1,col="white")
                    if (!is.null(folder)) {
                        filename <- sprintf("%s%04d.png", movie, movcount)
                        movcount <- movcount+1
                        rgl.snapshot(filename,fmt="png")
                    }
                }
            
            
            time1 <- Sys.time()
            if (!silent) {
                cat(paste("completed iteration",i, "in", round(time1-time0,2), "seconds\n"))
                cat("****************\n")
            }
        }
        invisible(mesh1)
    }


