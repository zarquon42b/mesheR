AmbergRegister <- function(mesh1,mesh2,lm1,lm2,k=1,lambda=1,iterations=15,rho=pi/2,dist=2,border=F,smooth=T)
    {
        mesh1rot <- rotmesh.onto(mesh1,lm1,lm2)
        lm1 <- mesh1rot$yrot
        mesh1 <- mesh1rot$mesh
        lmtmp1 <- lm1
        lmtmp2 <- lm2
        verts0 <- vert2points(mesh1)

        if (length(lambda) == 1)
            lambda <- rep(lambda,iterations)
        
        else if (length(lambda) != iterations)
            stop("lambda must be vector of length 'iterations'")

        tmp <- AmbergDeformSpam(mesh1,lmtmp1,lmtmp2,k0=k,lambda=lambda[1])
        if (iterations > 1)
            {
                for (i in 2:(iterations))
                    {
                        cat(paste("performing iteration",i,"\n"))
                        clost <- closemeshKD(tmp$mesh,mesh2)
                        verts1 <- vert2points(clost)
                        nc <- normcheck(clost,tmp$mesh)
                       

                        ## find valid hits
                        normgood <- as.logical(nc < rho)
                        distgood <- as.logical(abs(clost$quality) <= dist)
                        bordergood <- 1
                        if (border)
                            {
                                meshbord <- vcgBorder(mesh2)
                                bordgood <- as.logical(!meshbord$borderit[clost$ptr])
                            }
                        
                        good <- sort(which(as.logical(normgood*distgood*bordergood)))
                        lmtmp1 <- verts0[good,]
                        lmtmp2 <- verts1[good,]
                        tmp <- AmbergDeformSpam(mesh,lmtmp1,lmtmp2,k0=k,lambda=lambda[i],S=tmp$S)
                        if (smooth)
                            tmp$mesh <- vcgSmooth(tmp$mesh,iteration = 1)

                    }
            }
        return(tmp$mesh)
    }
