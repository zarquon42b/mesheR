gaussDisplace <- function(mesh1,mesh2,sigma,gamma=2,W0,f,oneway=F,k=1,nh=NULL,tol=0,cores=detectCores(),pro=c("morpho","vcg"),k0=50,prometh=1,rhotol=NULL,...)
{
### the workhorse function running in each iteration of gaussDisplMesh3d
  ## set projection function according to input request
  pro <- substring(pro[1],1L,1L)
  if (pro == "v")
    {project3d <- vcgClost
   }
  if (pro == "m")
    {
      protmp <- function(x,y,sign=F)
        {
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
  if (oneway)
    {
      M <- vert2points(mesh2)
    }
  else
    {
      Mpro <- project3d(mesh2,mesh1,sign=F)
      M <- vert2points(Mpro)
     }
  ## get neighbourhood for each point to minimize calculation time
  if (!is.null (nh))
    {
      clostIndW <- mcNNindex(S,W0,k=nh,cores=cores)
      clostIndP <- mcNNindex(M,W0,k=nh,cores=cores)
    }
  rt0 <- rep(0,dim(S)[1])
  rt1 <- rep(0,dim(M)[1])
  if (!is.null(rhotol))
    {
      rc <- rhotol
      rt0 <- normcheck(mesh1,Spro)
      if (!oneway)
        rt1 <- normcheck(mesh2,Mpro)
    }
  t3 <- Sys.time()
  D1 <- S-S0
  D2 <- M-M0
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
  iter <-floor(nx/cores)    
  for (i in 1:(cores-1))
    {
      mclist[[i]] <- list()
      mclist[[i]][[1]] <- W0[(1:iter)+((i-1)*iter),]
      mclist[[i]][[2]] <-c((1:iter)+((i-1)*iter))
    }
  mclist[[cores]] <-   list()
  mclist[[cores]][[1]] <- W0[-c(1:((cores-1)*iter)),]
  mclist[[cores]][[2]] <- c(1:dim(W0)[1])[-c(1:((cores-1)*iter))]
  ## define function to be run in parallel
  displacefun <- function(x,...)
    {
      tmp0 <- .Fortran("displace_mesh_gauss",x[[1]],nrow(x[[1]]),S0,nrow(S0),M,nrow(M),D1,D2,sigma,gamma,oneway,clostIndW[x[[2]],],nh,clostIndP[x[[2]],],tol=tol,rt0,rt1,rc,PACKAGE="Morpho")[[1]]
      return(tmp0)
    }
  tmp <- mclapply(mclist,displacefun,mc.cores=cores)
  
  for (i in 1:cores)
    {
      out <- rbind(out,tmp[[i]])
    }
  ## time <- system.time(out <- .Fortran("displace_mesh_gauss",W0,nrow(W0),S0,nrow(S0),M,nrow(M),D1,D2,sigma,gamma,oneway,clostIndW,nh,clostIndP,tol=tol))
  #print(time)
  addit <- W0+out
  return(list(addit=addit))
}

gaussDisplMesh3d <- function(mesh1,mesh2,iterations=10,smooth=NULL,smoothit=10,smoothtype=c("taubin","laplace","HClaplace"),sigma=20,gamma=2,f=1.2,oneway=F,tol=0,lm1=NULL,lm2=NULL,icp=FALSE,icpiter=3,uprange=0.95,rhotol=1,nh=NULL,toldist=0,patch=NULL,repro=FALSE,cores=detectCores(),pro=c("morpho","vcg"),k0=50,prometh=1,angtol=NULL,...)
  {
     if (is.null(nh))
      {
        nh=ceiling(150/meshres(mesh1))
        cat(paste("\nneighbourhood is set to",nh,"\n***************\n"))
      }
    ## set projection function according to input request
    pro <- substring(pro[1],1L,1L)
    if (pro == "v")
      {project3d <- vcgClost
     }
    if (pro == "m")
      {
        protmp <- function(x,y,sign=F)
          {
            out <- closemeshKD(x,y,k=k0,sign=sign)
            return(out)
          }
        project3d <- protmp
      }
    ## set up patch 
    rescue <- FALSE
    if (!is.null(patch))
      { ## append landmarks to meshes vertices
        mesh1$vb <- cbind(mesh1$vb,rbind(t(patch),1))
        colnames(mesh1$vb) <- NULL
        mdim <- dim(mesh1$vb)
        cols <- c((mdim[2]+1-dim(patch)[1]):mdim[2])
        if (substr(smoothtype,1L,1L) %in% c("h","H") )
          rescue <- TRUE
      }
    ## do icp matching
    if (icp)
      {
        cat("performing icp matching\n")
        mesh1 <- icp(mesh1,mesh2,lm1=lm1,lm2=lm2,uprange=uprange,rhotol=rhotol,iterations=icpiter,pro=pro,k=k0)
      }
    ## elastic matching starts
    cat("starting elastic matching\n****************\n")
    for (i in 1:iterations)
      {
        time0 <- Sys.time()
        if (!is.null(smooth))
          { if (i %% smooth == 0)
              {
                if(rescue)#keep patch from becoming NaN
                  {
                    tmppatch <- mesh1$vb[1:3,cols]
                  }
                cat("smoothing step\n")
                mesh1 <- vcgSmooth(mesh1,type=smoothtype,iteration=smoothit)
                cat("smoothing finished\n")
                if (rescue)
                  {
                    mesh1$vb[1:3,cols] <- tmppatch 
                  }
              }
          }
        ## call the workhorse doing the displacement
        tmp <- gaussDisplace(mesh1,mesh2,sigma=sigma,gamma=gamma,f=f,W0=vert2points(mesh1),nh=nh,k=i,tol=toldist,cores=cores,pro=pro,k0=k0,prometh=prometh,rhotol=angtol)
        mesh1$vb[1:3,] <- t(tmp$addit)
        ## project the patch back on the temporary surface
        if (!is.null(patch) && repro)
          {
            mesh1$vb[1:3,cols] <- project3d(t(mesh1$vb[1:3,cols]),mesh1)$vb[1:3,]
          }
        
        time1 <- Sys.time()
        cat(paste("completed iteration",i, "in", round(time1-time0,2), "seconds\n"))
        cat("****************\n")
      }
    if (!is.null(patch))
      invisible(list(mesh=mesh1,patch=vert2points(mesh1)[cols,]))
    else
      invisible(mesh1)
  }

    
