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
#' @details This is a quite simple attempt to remove an intersection after the alignment of two meshes. This is achieved by translating \code{reference} along the average difference vectors of all vertices
#' 
#' @export
removeIntersect <- function(reference,target,stepsize=0.2,maxit=100,tol=1,outside=TRUE,rot=FALSE) {
    mesh <- FALSE
    if (inherits(reference,"mesh3d")) {
        reference.old <- reference
        reference <- vert2points(reference)
        mesh <- TRUE
    }
    if (rgl.cur()==0)
        open3d()
    else
        rgl.clear()
    shade3d(target,col=3)
     ## if (icp)
    reference <- vert2points(icp(reference,target,iterations=1,silent=TRUE))
    points3d(reference,size=10)
    final <- FALSE
    count <- 0
    while (!final && count < maxit) {
        if (rot)
            reference <- rotinplace(vert2points(vcgClostKD(reference,target)),reference)
        dists <- vcgClostKD(reference,target)
        if (outside)
            below <- which(dists$quality > tol)
        else
            below <- which(dists$quality < tol)
       
        final <- !length(below)
        if (!final) {
            #print(range(dists$quality[below]))
            displacement <- diffs <- vert2points(dists)[below,]-reference[below,]
            if (is.matrix(diffs))
                displacement <- colMeans(diffs)
            displacement <- displacement/base::norm(displacement,"2")
            reference <- t(t(reference)+displacement*mean(abs(dists$quality)[below])*stepsize)
            count <- count+1
            rgl.pop()
            points3d(reference,size=10)
            ##Sys.sleep(0.2)
            print(count)
        }
    }
    if (mesh)
        reference <- updateVertices(reference.old,reference)
    return(reference)
}
           
rotinplace <- function(x,y) {
    rot <- rotonto(x,y)
    out <- t(t(rot$Y)+rot$transy)
    return(out)
}
