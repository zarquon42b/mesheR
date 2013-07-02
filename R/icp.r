#' Iterative closest point matching between two triangular meshes.
#' 
#' performs rigid body transformations (and scaling if requested) to map a
#' reference mesh onto a target mesh
#' 
#' the registration is done by minimising squared distances between reference
#' vertices and closest points on the target surface (a.k.a. Procrustes
#' registration)
#' 
#' @param mesh1 object of class "mesh3d". 
#' @param mesh2 object of class "mesh3d". 
#' @param iterations integer
#' @param scale lgoical: scale reference mesh to optimize fit.
#' @param lm1 optional: kx3 matrix containing reference points on mesh1. 
#' @param lm2 optional: kx3 matrix containing reference points on mesh2.
#' @param uprange quantile of distances between vertices of mesh1 and closest
#' points on mesh2. All hit points on mesh2 farther away than the specified
#' quantile are not used for matching.
#' @param maxdist maximum distance for closest points to be included for matching. Overrides uprange, if specified.
#' @param rhotol maximum allowed angle of which normals are allowed to differ
#' between reference points and hit target points. Hit points above this
#' threshold will be neglected.
#' @param k integer: number of closest triangles on target that will be
#' searched for hits. Only used when pro="morpho".
#' @param reflection logical: allow reflection.
#' @param pro character: algorithm for closest point search "morpho" calls
#' closemeshKD from the package Morpho, while vcg calls vcgClost from mesheR.
#' If the region of the targetmesh is much smaller than the region of the
#' reference, "vcg" can be really slow. Otherwise very fast. "morpho" is the
#' stable but somewhat slower algorithm.
#' 
#' @return Returns the rotated mesh1.
#' @author Stefan Schlager
#' @seealso \code{\link{rotmesh.onto}}, \code{\link{rotonto}}
#' @references Zhang Z. 1994. Iterative point matching for registration of
#' free-form curves and surfaces International Journal of Computer Vision
#' 13:119-152.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data(nose)
#' warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
#' rotnose <- icp(warpnose.long,shortnose.mesh,lm1=longnose.lm,lm2=shortnose.lm,rhotol=0.7,uprange=0.9)
#' shade3d(rotnose,col=2,alpha=0.7)
#' shade3d(shortnose.mesh,col=3,alpha=0.7)
#' 
#' @export icp
icp <- function(mesh1, mesh2, iterations=3,scale=T, lm1=NULL, lm2=NULL, uprange=0.9, maxdist=NULL, rhotol=pi, k=50, reflection=FALSE,pro=c("morpho","vcg"))
  {

    mesh1 <- adnormals(mesh1)
    mesh2 <- adnormals(mesh2)
    pro <- substring(pro[1],1L,1L)
    if (pro == "v")
    {project3d <- vcgClost
   }
  if (pro == "m")
    {
      tmpfun <- function(x,y,sign=F)
        {
          out <- closemeshKD(x,y,k=k,sign=sign)
          return(out)
        }
      project3d <- tmpfun
    }
    if (!is.null(lm1))## perform initial rough registration
      {
        mesh1 <- rotmesh.onto(mesh1,lm1,lm2,scale=scale,reflection=reflection)$mesh
      }
    
    for( i in 1:iterations)
      {
        cat("*")
        proMesh <- project3d(mesh1,mesh2,sign=F) ## project mesh1 onto mesh2
        x1 <- vert2points(mesh1)
        x2 <- vert2points(proMesh)
        dists <- abs(proMesh$quality)
        
        ## check if normals angles are below rhotol
        
        normchk <- normcheck(mesh1,proMesh)
        goodnorm <- which(normchk < rhotol)
        x1 <- x1[goodnorm,]
        x2 <- x2[goodnorm,]
        dists <- dists[goodnorm]
        
        ## check distances of remaining points and select points
        if (is.null(maxdist))
            qud <- quantile(dists,probs=uprange)
        else
            qud <- maxdist
        good <- which(dists < qud)
        mesh1 <- rotmesh.onto(mesh1,x1[good,],x2[good,],scale=scale)$mesh
      }
    cat("\n")
    return(mesh1)
  }
#' compare normal directions between two states of a mesh
#'
#' compare normal directions between two states of a mesh
#' @param mesh1 triangular mesh
#' @param mesh2 triangular mesh
#' @param circle logical: which method to use calculating the angel
#'
#' @return numeric vector containing angles between corresponding normals
#' @export normcheck

normcheck <- function(mesh1,mesh2,circle=TRUE)
  {
    x1norm <- mesh1$normals[1:3,]
    x2norm <- mesh2$normals[1:3,]
    ln <- dim(x2norm)[[2]]
    circle <- as.integer(circle)
    normcheck <- rep(0,dim(x2norm)[2])
    storage.mode(normcheck) <- "double"
    storage.mode(x1norm) <- "double"
    storage.mode(x2norm) <- "double"
    storage.mode(ln) <- "integer"
    storage.mode(circle) <- "integer"
    normcheck <- .Fortran("angcheck",x1norm,ln,x2norm,normcheck,circle)[[4]]
    return(normcheck)
  }
