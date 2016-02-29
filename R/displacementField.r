#' evaluate a displacement field using gaussian smoothed interpolation of a discrete displacement field
#'
#' evaluate a displacement field using gaussian smoothed interpolation of a given discrete displacement field
#' @param dispfield displacement field created using \link{\code{createDisplacementField}}
#' @param points matrix or mesh3d at which to evaluate the interpolated displacement field
#' @param sigma sigma controls the weight of the neighbourhood by defining the standard-deviation for the gaussian smoothing
#' @param gamma dampening factor (displacement vectors will be divided by \code{gamma}
#' @param k integer: number of k closest points to evaluate.
#' @param threads integer: number of threads to use for computing the interpolation.
#' @return returns an interpolated displacement field of class \code{displacement_field} at the positions of \code{points}.
#' @examples
#' require(Rvcg);require(Morpho)
#' data(dummyhead)
#' humoff <- meshOffset(dummyhead.mesh,offset=5)
#' dispfield <- createDisplacementField(dummyhead.mesh,humoff)
#'
#' highres <- vcgSubdivide(dummyhead.mesh)
#' ifield <- interpolateDisplacementField(dispfield,highres,threads=2,sigma = 50,k=500)
#' @export
interpolateDisplacementField <- function(dispfield, points, k=10, sigma=20,gamma=1, threads=parallel::detectCores()) {
    if (!inherits(dispfield,"displacement_field"))
        stop("please enter displacementfield")
    if (inherits(points,"mesh3d"))
        points <- vert2points(points)
    if (is.vector(points))
        points <- t(points)
    ## get closest points on domain
    clost <- vcgKDtree(dispfield$domain,points,k=k)
    clost$index <- clost$index-1L
    tmp = .Call("smoothField",points, dispfield$displacement_field,sigma,gamma,clost$index, clost$distance,threads)
    out <- list(domain=points,displacement_field=tmp)
    class(out) <- "displacement_field"
    return(out)
    
    
}
#' create a discrete displacement field
#'
#' create a discrete displacement field based on two sets of coordinates/meshes
#' @param reference k x 3 matrix containing coordinates or a mesh of classe "mesh3d"
#' @param target k x 3 matrix or a mesh of classe "mesh3d"
#' @return returns an object of class "displacement_field", which is a list containing
#' \item{domain}{k x matrix containing starting points of the displacement field}
#' \item{displacement_field}{k x matrix containing direction vectors of the displacement field}
#' @examples
#' require(Rvcg)
#' data(dummyhead)
#' humoff <- meshOffset(dummyhead.mesh,offset=5)
#' dispfield <- createDisplacementField(dummyhead.mesh,humoff)
#' @export
createDisplacementField <- function(reference,target) {
    if (inherits(reference,"mesh3d"))
        reference <- vert2points(reference)
     if (inherits(target,"mesh3d"))
         target <- vert2points(target)
    out <- list(domain=reference,displacement_field=target-reference)
    class(out) <- "displacement_field"
    return(out)
}

#' apply a discrete displacement field to a set of points/mesh in its domain
#'
#' apply a discrete displacement field to a set of points/mesh in its domain by applying the gaussian smoothed interpolation based of k closest neighbours
#' @param dispfield displacement field of class \code{displacement_field}
#' @param points matrix or mesh3d at which to evaluate the interpolated displacement field
#' @param sigma sigma controls the weight of the neighbourhood by defining the standard-deviation for the gaussian smoothing
#' @param gamma dampening factor (displacement vectors will be divided by \code{gamma}
#' @param k integer: number of k closest points to evaluate.
#' @param threads integer: number of threads to use for computing the interpolation.
#' @param interpolate logical:
#' @note if points is identical to the domain of the displacement field, no interpolation will be performed.
#' @return returns the displaced version of points
#' @export
applyDisplacementField <- function(dispfield,points,k=10,sigma=20,gamma=1, threads=1) {
    if (!inherits(dispfield,"displacement_field"))
         stop("please enter displacementfield")
    ## check if we need to interpolate at all
    if (!checkDispFieldDomain(dispfield,points))
        displacement <- interpolateDisplacementField(dispfield,points=points,k=k,sigma=sigma,gamma=gamma,threads=threads)
    else
        displacement <- dispfield
    updatePos <- displacement$domain+displacement$displacement_field
    if (inherits(points,"mesh3d")) {
        points$vb[1:3,] <- t(updatePos)
        if (!is.null(points$normals))
            points <- vcgUpdateNormals(points)
        return(points)
    } else {
        return(updatePos)
    }
   
    
}
#' visualize a displacement field
#'
#' visualize a displacement field
#' @param dispfield displacement field of class \code{displacement_field}
#' @param lwd width of the displacement vectors
#' @param colored logical: if TRUE, displacement vectors will be colored according to a heatmap.
#' @importFrom Morpho plotNormals meshDist
#' @importFrom colorRamps blue2green2red
#' @importFrom rgl wire3d
#' @export
plotDisplacementField <- function(dispfield,lwd=1,colored=TRUE) {
     if (!inherits(dispfield,"displacement_field"))
         stop("please enter displacementfield")
     if (!colored) {
     tmp <- list(vb=t(dispfield$domain),normals=t(dispfield$displacement_field))
     class(tmp) <- "mesh3d"
     plotNormals(tmp,lwd=lwd)
     } else {
         dists <- apply(dispfield$displacement_field,1,norm,"2")
         ramp <- colorRamps::blue2green2red(19)
         from <- 0
         to <- ceiling(max(dists))
         colseq <- seq(from=from,to=to,length.out=20)
         coldif <- colseq[2]-colseq[1]
         distqual <- ceiling((dists/coldif)+1e-14)
         colorall <- ramp[distqual]
         start <- dispfield$domain
         end <- dispfield$domain+dispfield$displacement_field
         vl <- nrow(start)
         dismesh <- list();class(dismesh) <- "mesh3d"
         dismesh$vb <- rbind(t(rbind(start,end)),1)
         dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
         dismesh$material$color <- rbind(colorall,colorall,colorall)
         wire3d(dismesh,lit=FALSE,lwd=lwd)
     }
}

## checks whether a set of points is identical to the domain of a displacement field
checkDispFieldDomain <- function(dispfield,points,tol=1e-12) {
    if (inherits(points,"mesh3d"))
        points <- vert2points(points)
    ndisp <- nrow(dispfield$domain)
    npoints <- nrow(points)
    if (ndisp != npoints)
        return(FALSE)
    else {
        chk <- vcgKDtree(dispfield$domain,points,k=1)
        if (isTRUE(all.equal(as.vector(chk$index),seq.int(ndisp),check.attributes = FALSE)))
            return(FALSE)
        if (max(chk$distance) > tol)
            return(FALSE)
    }
    return(TRUE)
}
