posteriorMeanICP <- function(model,target,rho=pi/2,subsample=500,dist=5,iterations=10,threads=0,test=T) {
    count <- 0
    ssmean <- initmesh <- vcgUpdateNormals(DrawMean(model))
    mysample <- 1:ncol(ssmean$vb)
    if (subsample > 0)
         mysample <- fastKmeans(initmesh,k=subsample,threads=threads)$selected
    while (count < iterations) {
        time0 <- Sys.time()
        vert_old <- vert2points(initmesh)
       
        useverts <- vert2points(initmesh)
        if (subsample > 0)
            useverts <- useverts[mysample,]
        clost <- vcgClostKD(useverts,target,k=50,threads=threads)
        verts1 <- vert2points(clost)
        if (subsample > 0) {
            cmesh <- initmesh
            cmesh$vb <- initmesh$vb[,mysample]
            cmesh$normals <- cmesh$normals[,mysample]
            nc <- normcheck(clost,cmesh)
            } else
                nc <- normcheck(clost,initmesh)                        
        
        ## find valid hits
        normgood <- as.logical(nc < rho)
        distgood <- as.logical(abs(clost$quality) <= dist)
        bordergood <- 1
        ## if (!border) 
        ##     bordgood <- as.logical(!meshbord$borderit[clost$faceptr])
                                        #dupes <- !(as.logical(vcgClean(clost)$remvert))
        dupes <- TRUE
        good <- sort(which(as.logical(normgood*distgood*bordergood*dupes)))
        if (test) {
            cmod <- statismoConstrainModel(model,pt=vert2points(ssmean)[mysample[good],],sample=verts1[mysample[good],],ptValueNoise=abs(clost$quality[mysample[good]]),computeScores = F)
        ## initmesh <- vcgUpdateNormals(PredictSample(model,lmModel=vert2points(ssmean)[good,],lmDataset=verts1[good,],align=FALSE,ptValueNoise=abs(clost$quality[good]),posteriorMean=T))
            initmesh <- vcgUpdateNormals(DrawMean(cmod))
        } else {
            initmesh <- vcgUpdateNormals(PredictSample(model,lmModel=vert2points(ssmean)[mysample[good],],lmDataset=verts1[good,],align=FALSE,ptValueNoise=1,posteriorMean=T))
        }
        time1 <- Sys.time()
        count <- count+1
        cat(paste("-> finished iteration",count,"in",round(as.numeric(time1-time0,unit="secs"),2), "seconds\n"))
    }
    return(initmesh)
}
        
