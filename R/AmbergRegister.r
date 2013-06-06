AmbergRegister <- function(mesh1,mesh2,lm1=NULL,lm2=NULL,k=1,lambda=1,iterations=15,rho=pi/2,dist=2,border=FALSE,smooth=TRUE,tol=1e-4,useiter=TRUE,minclost=50,distinc=1)
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
                mesh1rot <- rotmesh.onto(mesh1,lm1,lm2)
                lm1 <- mesh1rot$yrot
                meshorig <- mesh1 <- mesh1rot$mesh
                lmtmp1 <- lm1
                lmtmp2 <- lm2
                verts0 <- vert2points(mesh1)
                cat(paste("-> performing landmark based matching 1\n"))
                tmp <- AmbergDeformSpam(mesh1,lmtmp1,lmtmp2,k0=k[1],lambda=lambda[1])
                count <- count+1
                if (iterations == 1)
                    stopit <- TRUE
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
                        
                        clost <- closemeshKD(tmp$mesh,mesh2)
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
                            tmp$mesh <- vcgSmooth(tmp$mesh,iteration = 1)
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
