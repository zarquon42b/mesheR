icp <- function(mesh1,mesh2,iterations=3,scale=T,lm1=NULL,lm2=NULL,uprange=0.9,rhotol=pi,k=50,reflection=FALSE,pro=c("morpho","vcg"))
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
        qud <- quantile(dists,probs=uprange)
        good <- which(dists < qud)
        mesh1 <- rotmesh.onto(mesh1,x1[good,],x2[good,],scale=scale)$mesh
      }
    cat("\n")
    return(mesh1)
  }
        

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
    normcheck <- .Fortran("angcheck",x1norm,ln,x2norm,normcheck,circle,PACKAGE="Morpho")[[4]]
    return(normcheck)
  }
