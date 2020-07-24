#' This functions aims to translate two meshes to remove an intersection
#'
#' This functions aims to translate two meshes to remove an intersection. This is done based on signed distance
#'
#' @param reference matrix or mesh
#' @param target mesh of class mesh3d
#' @param stepsize regularize the actual displacement
#' @param maxit maximum iterations
#' @param tol positive number: stop if intersection is less than \code{tol}
#' @param outside logical: if TRUE the reference will be placed outside of the target. Inside otherwise
#' @param gradthresh numerical: (abs) gradient threshold to determine convergence
#' @param gradn integer: number of steps to include in computing gradient
#' @param realign logical: if TRUE, \code{reference} will be aligned by the points closest outside/inside of the target surface.
#' @param alignthresh threshold which points to consider for realignment. For \code{alignthresh=.5}, all coordinates within .5 mm outside/inside of the targtet will be used to realign the shape.
#' @param visualize logical: if TRUE watch the approximation
#' @param silent logical: suppress messages
#' @details This is a quite simple attempt to remove an intersection after the alignment of two meshes. This is achieved by translating \code{reference} along the average difference vectors of all vertices
#' @importFrom pracma gradient
#' @importFrom rgl rgl.cur rgl.clear
#' @export
removeIntersect <- function(reference,target,stepsize=0.2,maxit=100,tol=1,outside=TRUE,gradthresh=-Inf,gradn=Inf,realign=FALSE,alignthresh=.5,visualize=FALSE,silent=FALSE) {
    mesh <- FALSE
    ref.orig <- reference
    if (inherits(reference,"mesh3d")) {
         reference <- vert2points(reference)
         reference_old <- reference
         mesh <- TRUE
    }
    if (visualize) {
        if (rgl.cur()==0)
            open3d()
        else
            rgl.clear()
        shade3d(target,col=3)
        points3d(reference)
    }
    if (interactive() && visualize) 
        readline("please select viewpoint\n")
    final <- FALSE
    count <- 0
    tmpcount <- 1
    distvals <- NULL
    mygrad <- -Inf
    while (!final && count < maxit && abs(mygrad) > gradthresh) {
        dists <- vcgClostKD(reference,target)
        if (outside)
            below <- which(dists$quality > tol)
        else
            below <- which(dists$quality < tol)
        final <- !length(below)
        if (!silent) {
        if (length(below))
            cat(paste0("Max Dist:", max(abs(dists$quality[below])),"\n"))
        else
            cat(paste0("No points above the thresold of", tol," mm\n"))
        }
        if (!final) {
            distvals <- c(distvals,max(abs(dists$quality[below])))

            if (count > 1 && gradthresh > -Inf) {
                if (count > (gradn))
                    evaln <- (count-gradn+1):count+1
                else
                    evaln <- 1:length(distvals)

                f <- pracma::gradient(distvals[evaln])
                if (!silent)
                    cat(paste0("Gradient = ",f[length(f)],"\n"))
                mygrad <- f[length(f)]
            }
            
            displacement <- diffs <- vert2points(dists)[below,]-reference[below,]
            if (is.matrix(diffs))
                displacement <- colMeans(diffs)
            displacement <- displacement/base::norm(displacement,"2")
            reference <- t(t(reference)+displacement*mean(abs(dists$quality)[below])*stepsize)
            if (realign) {
                ## suppressMessages(reference <- vert2points(icp(reference,meshOffset(target,max(abs(dists$quality)[below])*inflate),iterations = 1,silent = TRUE)))
                myclost <- vcgClost(reference,target,sign=T)
                if (!outside)
                    myclost$quality <- -myclost$quality
                clost <- which(as.logical((myclost$quality < alignthresh)*(myclost$quality > 0)))
                
                #print(clost)
                if (length(clost)) {
                    myclostvert <- reference[clost,]
                    myclostaligned <- icp(myclostvert,target,iterations=1,getTransform = T)
                    reference <- applyTransform(reference,myclostaligned$transform)
                }
                
            }
            ## reference <- applyTransform(reference,myicp$transform)
            count <- count+1
            tmpcount <- tmpcount+1
            if (visualize) {
                rgl.pop()
                points3d(reference)
            }
            if (!silent)
            cat(paste0("finished iteration ",count,"\n\n"))
        }
    }
    if (mesh)
        reference <- updateVertices(ref.orig,reference)
    return(reference)
}

rotinplace <- function(x,y) {
    rot <- rotonto(x,y)
    out <- t(t(rot$Y)+rot$transy)
    return(out)
}







