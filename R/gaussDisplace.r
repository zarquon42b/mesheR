## @export gaussDisplace
gaussDisplace <- function(mesh1,mesh2,sigma,gamma=2,W0,f,oneway=F,k=1,nh=NULL,tol=0,cores=detectCores(),pro=c("morpho","vcg"),k0=50,prometh=1,rhotol=NULL,border=FALSE,horiz.disp=NULL,...)
{
### the workhorse function running in each iteration of gaussDisplMesh3d
    ## set projection function according to input request
    pro <- substring(pro[1],1L,1L)
    if (pro == "v") {
        project3d <- vcgClost
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
        clostIndW <- mcNNindex(S,W0,k=nh,cores=cores,...)
        if (!oneway)
            clostIndP <- mcNNindex(M,W0,k=nh,cores=cores,...)
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
    mclist <- list()
    nx <- dim(W0)[1]
    if (cores > 1) {
        iter <-floor(nx/cores)    
        for (i in 1:(cores-1)) {
            mclist[[i]] <- list()
            mclist[[i]][[1]] <- W0[(1:iter)+((i-1)*iter),]
            mclist[[i]][[2]] <-c((1:iter)+((i-1)*iter))
        }
        mclist[[cores]] <-   list()
        mclist[[cores]][[1]] <- W0[-c(1:((cores-1)*iter)),]
        mclist[[cores]][[2]] <- c(1:dim(W0)[1])[-c(1:((cores-1)*iter))]
    } else {
        mclist[[1]] <- list()
        mclist[[1]][[1]] <- W0
        mclist[[1]][[2]] <- 1:nx
    }
    ## define function to be run in parallel
    displacefun <- function(x,...)
        {
            tmp0 <- .Fortran("displace_mesh_gauss",x[[1]],nrow(x[[1]]),S0,nrow(S0),M,nrow(M),D1,D2,sigma,gamma,oneway,clostIndW[x[[2]],],nh,clostIndP[x[[2]],],tol=tol,rt0,rt1,rc,PACKAGE="mesheR")[[1]]
            return(tmp0)
        }
    tmp <- mclapply(mclist,displacefun,mc.cores=cores)
    
    for (i in 1:cores)
        out <- rbind(out,tmp[[i]])
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
#' @param mesh1 An object of class mesh3d used as atlas mesh to be deformed
#' onto the target mesh. Mesh resolution should be ~1.5.
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
#' @param icp Logical: if TRUE, iterative closest point procedure will be
#' executed.
#' @param icpiter Integer: Number of iterations of icp.
#' @param uprange argument passed to icp (see \code{\link{icp}})
#' @param rhotol Numeric: argument passed to \code{\link{icp}}.Exclude target
#' points with deviation of normals larger than than rhotol.
#' @param nh Integer: neighbourhood (number vertices) for controlling
#' displacement smoothing, default is 150/mesh resolution.
#' @param toldist Integer: Exclude everything from the whole procedure with a
#' greater distance from initial point than toldist. 0 disables this feature.
#' @param patch A m x 3 matrix containing the atlas landmark configuration on
#' mesh1 which is not present in mesh2 and will automatically placed on the
#' latter.
#' @param repro Logical: If TRUE, a reprojection of patch onto the iteratively
#' estimated target surface will be performed after each iteration.
#' @param cores integer: number of CPU cores to be used for multithreaded
#' subroutines.
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
#' @param Amberg vector containing 2 arguments invoking a smooth Deformation using \code{\link{AmbergDeformSpam}}. Layout: Amgerg=c(lambda, k0)
#' @param silent logical suppress messages
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
#' match <- gaussMatch(shortnose.mesh,warpnose.long,gamma=4,iterations=3,smooth=1,smoothtype="h",smoothit=10,nh=50,angtol=pi/2)
#' 
#' @export gaussMatch
#' @useDynLib mesheR
gaussMatch <- function(mesh1,mesh2,iterations=10,smooth=NULL,smoothit=10,smoothtype=c("taubin","laplace","HClaplace"),sigma=20,gamma=2,f=1.2,oneway=F,lm1=NULL,lm2=NULL,icp=FALSE,icpiter=3,uprange=0.95,rhotol=1,nh=NULL,toldist=0,patch=NULL,repro=FALSE,cores=detectCores(),pro=c("morpho","vcg"),k0=50,prometh=1,angtol=NULL,border=FALSE,horiz.disp=NULL,Amberg=NULL,silent=FALSE, ...)
    {
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
            project3d <- vcgClost
        } else if (pro == "m") {
            protmp <- function(x,y,sign=F) {
                out <- closemeshKD(x,y,k=k0,sign=sign)
                return(out)
            }
            project3d <- protmp
        }
        ## set up patch 
        rescue <- FALSE
        if (!is.null(patch)) { ## append landmarks to meshes vertices
            mesh1$vb <- cbind(mesh1$vb,rbind(t(patch),1))
            colnames(mesh1$vb) <- NULL
            mdim <- dim(mesh1$vb)
            cols <- c((mdim[2]+1-dim(patch)[1]):mdim[2])
            if (sum(substr(smoothtype,1L,1L) %in% c("h","H")) > 0 )
                rescue <- TRUE
        }
        ## do icp matching
        if (icp) {
            if (!silent)
                cat("performing icp matching\n")
            mesh1 <- icp(mesh1,mesh2,lm1=lm1,lm2=lm2,uprange=uprange,rhotol=rhotol,iterations=icpiter,pro=pro,k=k0)
        }
        ## elastic matching starts
        if (!silent)
            cat("starting elastic matching\n****************\n")
        for (i in 1:iterations) {
            time0 <- Sys.time()
            if (!is.null(smooth) && i > 1) {
                if (i %% smooth == 0) {
                    if(rescue && !is.null(patch))#keep patch from becoming NaN
                        tmppatch <- mesh1$vb[1:3,cols]
                    if (!silent)
                        cat("smoothing step\n")
                    mesh1 <- vcgSmooth(mesh1,type=smoothtype,iteration=smoothit)
                    if (!silent)
                        cat("smoothing finished\n")
                    if (rescue && !is.null(patch))
                        mesh1$vb[1:3,cols] <- tmppatch 
                }
            }
            ## call the workhorse doing the displacement
            tmp <- gaussDisplace(mesh1,mesh2,sigma=sigma,gamma=gamma,f=f,W0=vert2points(mesh1),nh=nh,k=i,tol=toldist,cores=cores,pro=pro,k0=k0,prometh=prometh,rhotol=angtol,border=border,oneway=oneway,horiz.disp = horiz.disp,...)
            
            if (!is.null(Amberg)) {#smooth deformation
                tmpmesh <- mesh1
                tmpmesh$vb[1:3,] <- t(tmp$addit)
                tmpmesh <- adnormals(mesh1)
                mesh1 <- AmbergDeformSpam(mesh1,vert2points(mesh1),tmp$addit,lambda=Amberg[1],k0=Amberg[2])$mesh
            } else
                mesh1$vb[1:3,] <- t(tmp$addit)
            mesh1 <- adnormals(mesh1)
            
            ## project the patch back on the temporary surface
            if (!is.null(patch) && repro)
                mesh1$vb[1:3,cols] <- project3d(t(mesh1$vb[1:3,cols]),mesh1)$vb[1:3,]
            
            time1 <- Sys.time()
            if (!silent) {
                cat(paste("completed iteration",i, "in", round(time1-time0,2), "seconds\n"))
                cat("****************\n")
            }
        }
        if (!is.null(patch))
            invisible(list(mesh=mesh1,patch=vert2points(mesh1)[cols,]))
        else
            invisible(mesh1)
    }


