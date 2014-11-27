outsideBBox <- function(x, corners) {
    box <- makeBox(corners)
    clost <- vcgClost(x,box,sign=T)
    inside <- which(clost$quality < 0)
    
}

makeBox <- function(corners) {
    bbox4 <- list(vb = rbind(t(corners),1))
    it <- c(1,2,4,4,3,1,1,5,2,2,5,6,7,8,6,6,5,7,7,3,4,4,8,7,1,3,7,7,5,1,2,6,4,4,6,8)
    
    bbox4$it <- matrix(it,3,length(it)/3)
    class(bbox4) <- "mesh3d"
    return(bbox4)
}
getMeshBox <- function(mesh,tri=FALSE,extend=0) {
    vv2 <- vert2points(mesh)
    pp <- prcomp(vv2)
    ranges <- apply(pp$x,2,range)
    ranges <- apply(ranges,2,extendrange,f=extend)
    mc <- as.matrix(expand.grid(ranges[,1],ranges[,2],ranges[,3]))
    bbox <- sweep(mc%*%t(pp$rotation),2,-pp$center)
    if (tri)
        bbox <- makeBox(mc)
    return(bbox)
}
