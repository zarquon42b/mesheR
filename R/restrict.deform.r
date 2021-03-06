#' restrict a landmark config or mesh to probabilistic bounds from a sample of landmarks/vertices
#'
#' restrict a landmark config or mesh to probabilistic bounds from a sample of landmarks/vertices
#' @param x matrix or triangular mesh of class "mesh3d"
#' @param model object of class "nosymproc" (output from \code{\link{procSym}} run on a set of corresponding landmarks/vertices), serving as a priori knowledge about vertex/landmark distribution.
#' @param sd standard deviation of PCscores defining boundaries
#' @param maxVar percentage to be explained by the PCs used for the restriction. E.g. if you want to include PCs explaining 95\% of the variance, set maxVar=95.
#' @param scale logical: if FALSE, configuration will be restricted within a probablistic hypercube, otherwise Chi-Square distribution will be used.
#' @param nPC number of PCs used in restricting (overrides maxVar)
#' @param probab if TRUE only the probability will be returned.
#' @param reference integer vector specifiying those landmarks that should be set to meanshape
#' @return
#' restricted landmarks/mesh
#'
#' @importFrom Morpho rotonto rotreverse vert2points cSize
#' @export
restrict <- function(x,model,sd=3,maxVar=95,scale=FALSE,nPC=NULL,probab=FALSE,reference=NULL) UseMethod("restrict")

#' @rdname restrict
#' @export
restrict.mesh3d <- function(x,model,sd=3,maxVar=95,scale=FALSE,nPC=NULL,probab=FALSE,reference=NULL) {
    mesh <- x
    x <- vert2points(mesh)
    out <- restrict(x,model=model,sd=sd,maxVar=maxVar,scale=scale,nPC=nPC,probab=FALSE,reference=NULL)
    mesh$vb[1:3,] <- t(out)
    mesh <- vcgUpdateNormals(mesh)
    return(mesh)
}
#' @rdname restrict
#' @export
restrict.matrix <- function(x,model,sd=3,maxVar=95,scale=FALSE,nPC=NULL,probab=FALSE,reference=NULL) {   
    dims <- dim(x)
    modAtt <- attributes(model)
    CSinit <- modAtt$CSinit
    rotscale <- modAtt$scale
    mshape <- model$mshape
    PCs <- model$PCs
    xorig <- x
    if (CSinit) {
        orsize <- cSize(x)
        x <- x/orsize
    }
    
    xrot <- rotonto(model$mshape,x,scale=rotscale)
    x <- xrot$yrot
    restr.x <- NULL
    if (is.null(nPC)) { ### select first # of PCs below threshold of total variance
        pc.used <- which(model$Variance[,3] < maxVar)
    } else {
        pc.used <- 1:nPC ## use predefined # of PCs
    }
    sds <-model$Variance[pc.used,1]
    prob <- TRUE
    sdl <- length(pc.used)
    xtmp <- x-mshape
    if (!is.null(reference)) {       
        xtmp[reference,] <- 0 ## set reference=mshape
    }
    xscore <- t(PCs[,pc.used])%*%as.vector(xtmp)
    
    if (scale) { ### use chisquare distribution of mahalanobis distance
        ##project into inverted eigenspace
        Mt <- qchisq(1-2*pnorm(sd,lower.tail=F),df=sdl)
        xscorePro <- xscore/sqrt(sds)
        probs <- sum(xscorePro^2)
        if (probs > Mt ) {
            prob=FALSE
            sca <- Mt/probs
            xscore <- xscorePro*sqrt(sca)*sqrt(sds)
        }
    } else { ### use probability hypercuboid
        sq.sds <- sqrt(sds)        
        for (i in 1:length(xscore)) {## check if PCscores are below threshold
            signum <- sign(xscore[i])
            if (abs(xscore[i]) > (sd*sq.sds[i])) {
                prob=FALSE
                xscore[i] <- sd*sq.sds[i]*signum
            }
        }
    }
    if (!probab) {
        restr.x <- matrix(PCs[,pc.used]%*%xscore,dims[1],dims[2])+mshape
        if (!is.null(reference)) {
            restr.x[reference,] <- x[reference,]
        }
    }
    restr.x <- rotreverse(restr.x,xrot)
    if (CSinit) {
        restr.x <- (restr.x/cSize(restr.x))*orsize
    }
    return(restr.x)
    ##return(list(restr.x=restr.x,prob=prob))
}

warpRestrict <- function(x,which,tar.lm,model,tol=1e-5,sd=3,maxVar=95,scale=F,recurse=T,uniform=TRUE,iterations=NULL,nPC=NULL,stop.prob=TRUE,spline=TRUE,useReference=FALSE)
  {

    reference <- NULL
    if (useReference)### "which" will be ignored when calculating probability
      {
        reference <- which
      }
    tmp.res <- list()
    tmp.res$prob <-  FALSE
    if (is.null(nPC))
      {
        pc.used <-which(model$Variance[,3] < maxVar)
        print(paste("First ",max(pc.used)," PCs used"))
      }
    
    x.lm <- x[which,]
    sd.i <- sd    
    tmp <- x
    tmp.lm <- x.lm
    p <- 1e10
    count <- 0
    tmp.orig <- tmp
    tmp.old <- tmp     
    tmp <- tps3d(tmp,tmp.lm,tar.lm)## warp onto target
    tmp <- rotonto(model$mshape,tmp,scale=T)$yrot ### register in database space
    
    while(p > tol)
      {
        
        cat(paste("running iteration",count,"\n"))
        p.old <- p
        tmp.old <- tmp
        prob <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale,nPC=nPC,probab=T)$prob
        if (!prob) ### not yet probable
          {
            ## restrict to boundaries
           
            tmp <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale,nPC=nPC,probab=F,reference=reference)$restr.x
           # print(dim(tmp))
             tmp.orig <- tmp
            if (spline) ###use restricted data and warp it onto target
              {
                tmp.lm <- tmp[which,] 
                tmp <- tps3d(tmp,tmp.lm,tar.lm)
                tmp <- rotonto(model$mshape,tmp,scale=T)$yrot ### register in database space
              }
            else ### replace restricted reference lm with original ones
              {
                tmp[which,] <- tmp.orig[which,]
                tmp <- rotonto(model$mshape,tmp,scale=T)$yrot ### register in database space
              }
          }
        
        if (prob && stop.prob)
          {
            cat(paste("probable shape within the boundaries of",sd, "sd reached\n"))
                                        #if (!sd.i >= sd)
            p <- 0
          }
        else
          {
            p <- angle.calc(tmp,tmp.old)
          }                         # print(p)
        if (uniform)
          {
            if (p > p.old)
              {
                p <- 0
                tmp <- tmp.old
              }                  
          }
        if (!is.null(iterations))
          if (count == iterations)
            {
              p <- 0
            }
        count <- count+1
      }
    if (prob)
      {
        cat(paste("a probable shape was reached after",count-1,"iterations\n"))
      }
    else
      {
        cat("progress terminated without reaching probability\n")
      }
    tmp.orig <- rotonmat(tmp.orig,tmp.orig[which,],tar.lm,scale=T)
    tmp <- rotonmat(tmp,tmp[which,],tar.lm,scale=T)
    clean.out <- tmp
    
    return(list(raw=tmp.orig,clean=clean.out,tmp=tmp.res))
  }
