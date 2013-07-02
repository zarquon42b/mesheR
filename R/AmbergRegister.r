#' Register two triangular meshes based on smooth deformation.
#' 
#' Perform registration of two triangular meshes, minimizing per-face
#' distortions. 
#' 
#' @param mesh1 reference mesh: triangular mesh of class "mesh3d". No loose
#' vertices, edges and degenerated faces are allowed. 
#' @param mesh2 target mesh: triangular mesh of class "mesh3d". 
#' @param lm1 m x 3 matrix containing correspondences on "mesh1" 
#' @param lm2 m x 3 matrix containing target correspondences on "mesh2" 
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
#' @param scale logical: if TRUE, initial landmark based rigid registration
#' includes scaling. 
#' @param icp vector of length 4. Passing parameters to \code{\link{icp}},
#' which is performed after intial landmark based registration. The parameters
#' are icp[1]=iterations; icp[2]=rhotol; icp[3]=uprange, and icp[4]=scale. If
#' icp=NULL, no ICP-matching is performed.  E.g. icp=c(3,pi/2,0.6,TRUE) will
#' result in 3 icp iterations, condidering the closest 60\% of correspondences
#' with normal deviation of pi/2 and include scaling.
#' @param cores integer: how many cores to use for closest point search
#' @return 
#' \item{mesh}{registered mesh}
#' \item{meshrot }{mesh1, rotated onto mesh2}
#' \item{lm1rot }{lm1, rotated onto lm2}
#' \item{lmtmp1 }{correspondences on updated reference mesh of last iteration}
#' \item{lmtmp2 }{correspondences on updated target mesh of last iteration}
#' @author Stefan Schlager
#' @seealso \code{\link{gaussMatch}}
#' @references Amberg, B. 2011. Editing faces in videos, University of Basel.
#' @keywords ~kwd1 ~kwd2
#' @export AmbergRegister
AmbergRegister <- function(mesh1, mesh2, lm1=NULL, lm2=NULL, k=1, lambda=1, iterations=15, rho=pi/2, dist=2, border=FALSE, smooth=TRUE, smoothit=1, smoothtype="t", tol=1e-4, useiter=TRUE, minclost=50, distinc=1, scale=TRUE, icp=NULL, cores=1)
    {
        mesh1 <- rmUnrefVertex(mesh1)
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
        
        meshorig <- mesh1
        stopit <- FALSE
        if (!is.null(lm1) && !is.null(lm2))
            {   ## case: landmarks are provided
                if (!is.null(icp))##perform initial icp-matching
                    {
                        meshorig <- mesh1 <- icp(mesh1,mesh2,lm1=lm1,lm2=lm2,iterations=icp[1],rhotol=icp[2],uprange=icp[3],scale=icp[4])
                        tmp <- list()
                        tmp$mesh <- mesh1
                        if (!useiter)
                            tmp$S <- createS(mesh1)
                        verts0 <- vert2points(mesh1)
                                            }
                else
                    {
                        mesh1rot <- rotmesh.onto(mesh1,lm1,lm2,scale=scale)
                        lm1 <- mesh1rot$yrot
                        meshorig <- mesh1 <- mesh1rot$mesh
                        lmtmp1 <- lm1
                        lmtmp2 <- lm2
                        cat(paste("-> performing landmark based matching 1\n"))
                        tmp <- AmbergDeformSpam(mesh1,lmtmp1,lmtmp2,k0=k[1],lambda=lambda[1])
                        count <- count+1
                        if (iterations == 1)
                            stopit <- TRUE
                    }
                verts0 <- vert2points(mesh1)
            }
        else
            {   ## case: meshes are already aligned
                tmp <- list()
                tmp$mesh <- mesh1
                if (!useiter)
                    tmp$S <- createS(mesh1)
                verts0 <- vert2points(mesh1)
            }
        
        if (!stopit)
            {
                ## set error and counter appropriately
                error <- 1e12
                count <- count+1
                while (count <= iterations && error > tol)
                    {
                        time0 <- Sys.time()
                        if (useiter)
                            {
                                verts0 <- vert2points(tmp$mesh)
                                mesh1 <- tmp$mesh
                            }
                        vert_old <- vert2points(tmp$mesh)
                        
                        clost <- closemeshKD(tmp$mesh,mesh2, cores=cores)
                        verts1 <- vert2points(clost)

                        nc <- normcheck(clost,tmp$mesh)                        

                        ## find valid hits
                        normgood <- as.logical(nc < rho)
                        distgood <- as.logical(abs(clost$quality) <= dist)
                        bordergood <- 1
                        if (!border)
                            {
                                bordgood <- as.logical(!meshbord$borderit[clost$ptr])
                            }
                        
                        good <- sort(which(as.logical(normgood*distgood*bordergood)))
                        
### in case no good hit is found within the given distance we increase the distance by 1mm until valid references are found:
                        increase <- distinc
                        while (length(good) < minclost)
                            {
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
                        tmp <- AmbergDeformSpam(mesh1,lmtmp1,lmtmp2,k0=k[count],lambda=lambda[count],S=tmp$S)
                                        #oo <- wire3d(tmp$mesh,col=count)
                        gc()
                        ## calculate error
                        if (smooth)
                            tmp$mesh <- vcgSmooth(tmp$mesh,iteration = smoothit,type=smoothtype)
                        error <- sum((vert2points(tmp$mesh)-vert_old)^2)/nrow(vert_old)
                        time1 <- Sys.time()
                        cat(paste("-> finished iteration",count,"in",round(time1-time0,2), "seconds\n"))
                        cat(paste(" Info: MSE between iterations:",error,"\n"))
                        if (error < tol)
                            cat(paste("***\n==> Convergence threshold reached after",count,"iterations\n"))
                        count <- count+1
                    }
            }
        return(list(mesh=tmp$mesh,meshrot=meshorig,lm1rot=lm1,lmtmp1=lmtmp1,lmtmp2=lmtmp2))
    }
