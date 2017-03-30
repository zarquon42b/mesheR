#' Register two triangular meshes based on smooth deformation.
#' 
#' Perform registration of two triangular meshes, minimizing per-face
#' distortions. 
#' 
#' @param x reference mesh: triangular mesh of class "mesh3d"or of class BayesDeform created by createBayes to restrict based on a known distribution. To use this option the package RvtkStatismo \url{https://github.com/zarquon42b/RvtkStatismo} has to be installed.No loose vertices, edges and degenerated faces are allowed. 
#' @param mesh2 target mesh: triangular mesh of class "mesh3d". 
#' @param lm1 m x 3 matrix containing correspondences on "mesh1". 
#' @param lm2 m x 3 matrix containing target correspondences on "mesh2".
#' @param k integer: parameter regularizing face normal distortion. Can be
#' vector of length(iterations) or single value. 
#' @param lambda numeric: parameter regularizing faces's distortion. Can be
#' vector of length(iterations) or single value.
#' @param iterations integer: number of iterations to run. 
#' @param rho numeric: 0 < rho < 2*pi tolerance of normal deviation between
#' reference vertices and corresponding closest points on target suface. 
#' @param dist numeric: tolerance of maximal distance between reference
#' vertices and corresponding closest points on target suface.
#' @param border logical: if FALSE, hits on border faces are ignored (reduces
#' distortion) 
#' @param smooth logical: if TRUE after each iteration a mesh smoothing is performed.
#' @param smoothit integer: determine smoothing iterations.
#' @param smoothtype character: select smoothing algorithm - see vcgSmooth for further details.
#' @param tol numeric: convergence threshold of MSE between vertices of two
#' successive iterations.
#' @param useiter logical: if TRUE, each iteration uses the updated reference
#' mesh, if false. The original mesh will be deformed based on the updated
#' correspondences. 
#' @param minclost minimum amount of correspondence points. If less
#' correspondences are found, dist will be increased by "distinc" (see below).
#' @param distinc increment of dist, in case minclost is not reached.
#' @param rigid named list. Passing parameters to \code{\link{icp}}, for rigid registration. If landmarks are provided and only those should count, set rigid$iterations=0.
#' @param similarity named list. Passing parameters to \code{\link{icp}}, for similarity registration (rigid +scaling). If landmarks are provided and only those should count, set similarity$iterations=0 (and rigid=NULL).
#'@param affine named list. Passing parameters to \code{\link{icp}}, for affine registration. If landmarks are provided and only those should count, set similarity$iterations=0 (with rigid=NULL and similarity=NULL)
#'  @param tps logical: if TRUE and landmarks are provided, the reference will be mapped to the target using a Thin-Plate Spline interpolation. Overrides \code{rigid},\code{affine} and \code{similarity}.
#' @param pcAlign if TRUE, surfaces are prealigned by principal axis. Overrides intial landmark based alignment.
#' @param nn integer: closest barycenters. During search for closest points on target, the closest \code{nn} faces are probed. The larger \code{nn} is , the more accurate the closest point search but also the more time consuming. If landmarks are provided and only those should count, set rigid$iterations=0.
#' @param silent logical: no verbosity
#' @param useConstrained logical: if TRUE and Bayes and landmarks are defined, the landmarks are not only used to get a suitable reference but the model will also be constrained by the landmarks to subsequently restrict the shape variability. If FALSE, the full model is used.
#' @param forceLM logical: if icp is requested landmark based deformation will be applied after icp-based transformation.
#' @param visualize logical request visualization of deformation process.
#' @param folder logical: if visualize=TRUE, this can specify a folder to save screenshots of each deformation state, in order to create a movie or an animated gif.
#' @param noinc logical: if TRUE and x is of class 'Bayes', the process stops if the distance from the target to the deformed reference increases compared to the previous iteration.
#' @param bboxCrop extend of the bounding box around mesh1 (after alignmend) that will be cropped from target to speed things up.
#' @param threads integer: threads to use in closest point search.
#' @return 
#' \item{mesh}{registered mesh}
#' \item{affine }{affine 4x4 transformation matrix mapping mesh1 onto mesh2}
#' \item{lm1 }{lm1 mapped onto the registered template}
#' 
#' @details This function runs an elastic-ICP surface matching algorithm, that minimizes the original meshes internal structure by solving a sparse equation system. The user can control 2 parameters of mesh stiffness: \code{lambda} and \code{k}. \code{lambda} controls the impact of the control points (closest points) as it is a weight applied to the equation system. The value of \code{lambda} should be carefully selected depending on the object overall size: i.e. to match two tiny meshes one will need a higher value than a for a larger object (example: I found values between 0 and 1 suitable for human faces and values between 10 and 100 suitable for mice teeth). \code{k} controls the normal slackness, i.e. the deviation of normal direction. The larger, \code{k}, the more elastic the deformation will be. \code{lambda} and \code{k} can be specified as vectors of length \code{iterations}, to assign a specific value for each iteration. 
#' @author Stefan Schlager
#' @seealso \code{\link{gaussMatch}}
#' @references Amberg, B. 2011. Editing faces in videos, University of Basel.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' require(Morpho)
#' 
#' require(Rvcg)
#' data(humface)
#' data(dummyhead)
#' ## set parameters making each iteration more elastic
#' # only 10 iterations to keep example calculation time reasonable.
#' params <- list(iterations=10) 
#' params <- append(params, list(
#'    # first lambda is set relatively high because first matching uses landmarks
#'    # then let it increase from 0.2 to 0.6
#'    lambda=c(0.7,seq(from = 0.2,to=0.6,length.out = params$iterations-1)),
#'    # treat k similar as lambda
#'    k=c(10,seq(from = 1,to=params$iterations-1,by=1)),
#'    useiter=FALSE # iteratively deform dummyhead onto humface
#'    ))
#' #we also want the landmarks to be used in an initial similarity transform
#' similarity <- list(iterations=0)
#' map <- AmbergRegister(dummyhead.mesh, humface, lm1=dummyhead.lm,
#'                  lm2=humface.lm, iterations=params$iterations,similarity=similarity,
#'                  k=params$k, lambda=params$lambda, useiter=params$useiter)
#' # compare matched and original face:
#' \dontrun{
#' require(rgl)
#' meshDist(map$mesh, humface ,from=-3,to=3,tol=0.5)
#' # render original mesh as wireframe
#' shade3d(humface,front="lines",back="lines")
#' }
#' ##example with different icp matchings:
#' rigid <- list(iterations=30,subsample=200,rhotol=pi/2,uprange=0.6)
#' similarity <- list(iterations=30, subsample=200,rhotol=pi/2,uprange=0.6)
#' affine <- list(iterations=30,subsample=200,rhotol=pi/2,uprange=0.6)
#' map <- AmbergRegister(dummyhead.mesh, humface, lm1=dummyhead.lm,
#'                       lm2=humface.lm, iterations=params$iterations,
#'                       k=params$k, lambda=params$lambda, useiter=params$useiter,rigid=rigid,
#'                       similarity=similarity,affine=affine,forceLM = TRUE)
#' @importFrom Rvcg vcgClean vcgClost vcgUpdateNormals
#' @importFrom Morpho meshcube applyTransform computeTransform pcAlign
#' @export AmbergRegister
AmbergRegister <- function(x, mesh2, lm1=NULL, lm2=NULL, k=1, lambda=1, iterations=15, rho=pi/2, dist=2, border=FALSE, smooth=TRUE, smoothit=1, smoothtype="t", tol=1e-10, useiter=TRUE, minclost=50, distinc=1, rigid=NULL,similarity=NULL, affine=NULL,tps=FALSE, pcAlign=FALSE,nn=20, silent=FALSE, useConstrained=TRUE, forceLM=FALSE,visualize=FALSE, folder=NULL,noinc=FALSE,bboxCrop=NULL,threads=0)
{
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
    mesh1 <- rmUnrefVertex(mesh1, silent=TRUE)
    mesh1 <- vcgUpdateNormals(mesh1)
    mesh2 <- vcgUpdateNormals(mesh2)
    meshbord <- vcgBorder(mesh2)
    count <- 0
    if (iterations < 1)
        iterations <- 1e10
    if (length(lambda) == 1)
        lambda <- rep(lambda,iterations)
    else if (length(lambda) != iterations)
        stop("lambda must be vector of length 'iterations'")
    k <- round(k)# make sure k is integer - otherwise RAM overkill
    if (length(k) == 1)
        k <- rep(k,iterations)
    else if (length(k) != iterations)
        stop("k must be vector of length 'iterations'")

    affinemat <- NULL
    meshorig <- mesh1
    stopit <- FALSE
    hasLM <- FALSE
    if (!is.null(lm1) && !is.null(lm2)) {
        hasLM <- TRUE
        bary <- vcgClost(lm1,mesh1,barycentric = T)
    }
    if (!is.null(Bayes$initparams)) {
        mesh1 <- RvtkStatismo::DrawSample(Bayes$model,Bayes$initparams)
        if (hasLM)
            lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
                                        #wire3d(mesh1);spheres3d(lm1)
                                        #return(1)
    }
    if (hasLM) {## case: landmarks are provided
        
        if (!is.null(Bayes) && hasLM) {
            ##register landmarks on model and constrain reference
            lm2tmp <- rotonto(lm1,lm2,scale=Bayes$model@scale,reflection=FALSE)$yrot
            constMod <- RvtkStatismo::statismoConstrainModel(Bayes$model,lm2tmp,lm1,Bayes$ptValueNoise)
            if (useConstrained) {
                Bayes$model <- constMod
                mesh1 <- vcgUpdateNormals(RvtkStatismo::DrawMean(Bayes$model))
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
        }                
        if (tps) {
            mesh1 <- tps3d(mesh1,lm1,lm2,threads=threads)
            lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
        } else {
            if (is.null(rigid) && is.null(affine) && is.null(similarity) && !pcAlign) {
                if (is.null(Bayes)) {                
                    rigid <- list(iterations=0)
                    if (!silent)
                        cat("\n landmarks but no transform specified, performing rigid transform\n")
                } else if (Bayes$align) {
                    rigid <- list(iterations=0)
                    if (!silent)
                        cat("\n landmarks but no transform specified, performing rigid transform\n")
                }
            }
            if (pcAlign) {
                mesh1 <- pcAlign(mesh1,mesh2)
                lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            }
            if (!is.null(rigid)) { ##perform rigid icp-matching
                if (!pcAlign) {
                    rigid$lm1 <- lm1
                    rigid$lm2 <- lm2
                }
                mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent,threads=threads)
                if (hasLM)
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
        
        affinemat <- computeTransform(vert2points(mesh1),vert2points(meshorig))
        tmp <- list()
        tmp$mesh <- mesh1
        if (!useiter && !forceLM)
            tmp$S <- createS(mesh1)
        
        if (forceLM && hasLM) {
            lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
            tmp <- AmbergDeformSpam(mesh1,lm1,lm2,k0=k[1],lambda=lambda[1])
            count <- count+1
            if (iterations == 1)
                stopit <- TRUE
        }
        verts0 <- vert2points(mesh1)
        
        
    } else if (pcAlign || !is.null(rigid) || !is.null(affine) || !is.null(similarity)) {
        if (pcAlign) {
            mesh1 <- pcAlign(mesh1,mesh2)
        }
        if (!is.null(rigid)){ ##perform rigid icp-matching
            mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent)
        }
        if (!is.null(similarity)) {##similarity matching
            mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent)
        }
        if (!is.null(affine)) {##similarity matching
            mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent)
        }
        affinemat <- computeTransform(vert2points(mesh1),vert2points(meshorig))
        tmp <- list(mesh=mesh1)
        if (!useiter)
            tmp$S <- createS(mesh1)
        verts0 <- vert2points(mesh1)
        
    } else {   ## case: meshes are already aligned
        affinemat <- diag(4)
        tmp <- list()
        tmp$mesh <- mesh1
        if (!useiter)
            tmp$S <- createS(mesh1)
        verts0 <- vert2points(mesh1)
    }
    if (!is.null(bboxCrop)) {
        mesh2 <- cropOutsideBBox(mesh1,mesh2,extend=bboxCrop)
        if (!silent)
            cat("cropping target mesh\n")
    } 
    if (visualize) {
        rglid <- NULL
        if (!length(rgl.ids()$id)) 
            open3d()
        else {
            rgl.bringtotop()
            rgl.clear()
        }
        points3d(meshcube(tmp$mesh),col="white",alpha=0)
        shade3d(mesh2,col=2,specular=1)
        if (!is.null(rglid))
            rgl.pop(id=rglid)
        rglid <- shade3d(tmp$mesh,col="white",front="lines", back="lines")
        
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
    
    if (!stopit) {
        ## set error and counter appropriately
        distance <- 1e12
        error <- 1e12
        count <- count+1
        while (count <= iterations && error > tol) {
            time0 <- Sys.time()
            if (useiter) {
                verts0 <- vert2points(tmp$mesh)
                mesh1 <- tmp$mesh
            }
            vert_old <- vert2points(tmp$mesh)
            clost <- vcgClostKD(tmp$mesh,mesh2,k=nn,threads=threads)
            verts1 <- vert2points(clost)
            nc <- normcheck(clost,tmp$mesh,threads = threads)                        
            
            ## find valid hits
            normgood <- as.logical(nc < rho)
            distgood <- as.logical(abs(clost$quality) <= dist)
            bordergood <- 1
            if (!border) 
                bordgood <- as.logical(!meshbord$borderit[clost$faceptr])
                                        #dupes <- !(as.logical(vcgClean(clost)$remvert))
            dupes <- TRUE
            good <- sort(which(as.logical(normgood*distgood*bordergood*dupes)))
            
            
            
            
### in case no good hit is found within the given distance we increase the distance by 1mm until valid references are found:
            increase <- distinc
            while (length(good) < minclost) {
                distgood <- as.logical(abs(clost$quality) <= (dist+increase))
                good <- sort(which(as.logical(normgood*distgood*bordergood)))
                increase <- increase+distinc
                cat(paste("distance increased to",dist+increase,"\n"))
            }
            
            
            ## update reference points
            lmtmp1 <- verts0[good,]
            lmtmp2 <- verts1[good,]
            ## map it according to new reference points
                                        #points3d(lmtmp2,col=count)
            if (useiter)
                tmp$S <- NULL

            tmpold <- tmp
            chk <- try(tmp <- AmbergDeformSpam(mesh1,lmtmp1,lmtmp2,k0=k[count],lambda=lambda[count],S=tmp$S),silent = TRUE)
            if (inherits(chk,"try-error")) {
                tmp <- tmpold
                cat("iteration failed: previous iteration used")
            }
            gc()
            if (smooth)
                tmp$mesh <- vcgSmooth(tmp$mesh,iteration = smoothit,type=smoothtype)
            ## calculate error
            if (!is.null(Bayes) && length(Bayes$sdmax) >= count) {
                if (!is.null(Bayes$wt)) {
                    wt <- Bayes$wt[count]
                    wts <- c(1,wt)
                    wts <- wts/sum(wts)
                    tmpmesh <- RvtkStatismo::PredictSample(Bayes$model,tmp$mesh,TRUE, sdmax=Bayes$sdmax[count],align=Bayes$align,mahaprob=Bayes$mahaprob)
                    tmp$mesh$vb[1:3,] <- wts[1]*tmp$mesh$vb[1:3,]+wts[2]*tmpmesh$vb[1:3,]
                } else {
                    tmp$mesh <- RvtkStatismo::PredictSample(Bayes$model,tmp$mesh,TRUE, sdmax=Bayes$sdmax[count],align=Bayes$align,mahaprob=Bayes$mahaprob)
                }
                
            }
            distance_old <- distance
            distance <- mean(vcgClostKD(mesh2,tmp$mesh,k0=10,sign=F,threads=threads)$quality)
            if (distance > distance_old && !is.null(Bayes) && noinc) {
                cat("\n=========================================\n")
                message(paste(" Info: Distance is increasing, matching stopped after ",count,"iterations\n"))
                count <- 1e10
                tmp <- tmpold
            }
            tmp$mesh <- vcgUpdateNormals(tmp$mesh)
            if (visualize) {
                
                if (!is.null(rglid))
                    rgl.pop(id=rglid)
                rglid <- shade3d(tmp$mesh,col="white",front="lines",back="lines")
                if (!is.null(folder)) {
                    filename <- sprintf("%s%04d.png", movie, movcount)
                    movcount <- movcount+1
                    rgl.snapshot(filename,fmt="png")
                }
            }
            error <- sum((vert2points(tmp$mesh)-vert_old)^2)/nrow(vert_old)
            
            time1 <- Sys.time()
            if (!silent && count < 1e10) {
                cat(paste("-> finished iteration",count,"in",round(as.numeric(time1-time0,unit="secs"),2), "seconds\n"))
                cat(paste(" Info: MSE between iterations:",error,"\n"))
                cat(paste(" Info: Average distance to target:",distance,"\n"))
                if (error < tol)
                    cat(paste("***\n==> Convergence threshold reached after",count,"iterations\n"))
            }
            count <- count+1
        }
    }
    lm1map <- NULL
    if (!is.null(lm1))
        lm1map <- lm1 <- bary2point(bary$barycoords,bary$faceptr,tmp$mesh)
    return(list(mesh=tmp$mesh,affine=affinemat,lm1=lm1map))
}


rigSimAff <- function(mesh1,mesh2,args,type="r",silent=TRUE,threads=1) {
    iterations <- args$iterations; if (is.null(iterations)) iterations <- 3
    lm1=args$lm1
    lm2 <- args$lm2
    uprange <- args$uprange; if (is.null(uprange)) uprange <- 0.9
    maxdist <- args$maxdist
    minclost <- args$minclost; if (is.null(minclost)) minclost <- 50
    distinc <- args$distinc;
    rhotol <- args$rhotol; if (is.null(rhotol)) rhotol <- pi
    k <- args$k; if (is.null(k)) k <- 50
    reflection <- args$reflection;  if (is.null(reflection)) reflection <- FALSE
    subsample <- args$subsample
    out <- icp(mesh1, mesh2, iterations=iterations, lm1=lm1, lm2=lm2, uprange=uprange ,maxdist=maxdist, minclost=minclost, distinc=distinc, rhotol=rhotol, k=k, reflection=reflection, silent = silent,subsample=subsample, type=type,threads=threads)
    return(out)
}
