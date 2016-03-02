#' create a discrete displacement field
#'
#' create a discrete displacement field based on two sets of coordinates/meshes
#' @param reference k x 3 matrix containing coordinates or a mesh of classe "mesh3d"
#' @param target k x 3 matrix or a mesh of classe "mesh3d"
#' @return returns an object of class "DisplacementField", which is a list containing
#' \item{domain}{k x matrix containing starting points of the displacement field}
#' \item{DisplacementField}{k x matrix containing direction vectors of the displacement field}
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
    out <- list(domain=reference,DisplacementField=target-reference)
    class(out) <- "DisplacementField"
    return(out)
}

#' invert a displacement field
#'
#' invert a displacement field
#' @param dispfield displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @return returns an inverted displacement field
#' @examples
#' require(Rvcg);require(Morpho)
#' data(dummyhead)
#' humoff <- meshOffset(dummyhead.mesh,offset=5)
#' dispfield <- createDisplacementField(dummyhead.mesh,humoff)
#' dispfield_inverted <- invertDisplacementField(dispfield)
#' @export
invertDisplacementField <- function(dispfield) {
    validDisplaceField(dispfield)
    dispfield$domain <- dispfield$domain+dispfield$DisplacementField
    dispfield$DisplacementField <- -dispfield$DisplacementField
    return(dispfield)
}
    

#' evaluate a displacement field using gaussian smoothed interpolation of a discrete displacement field
#'
#' evaluate a displacement field using gaussian smoothed interpolation of a given discrete displacement field
#' @param dispfield @param dispfield displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @param k integer: number of k closest points to evaluate.
#' @param points matrix or mesh3d at which to evaluate the interpolated displacement field
#' @param sigma kernel bandwidth used for smoothing
#' @param gamma dampening factor (displacement vectors will be divided by \code{gamma}
#' @param type kernel function for smoothing are "Gauss","Laplace" and "Exponential"
#' @param threads integer: number of threads to use for computing the interpolation.
#' @return returns an interpolated displacement field of class \code{displacement_field} at the positions of \code{points}.
#' @note The k-closest coordinates of the displacement field are used to calculate a weighted (smoothed) displacement field for each point. The displacement field can then optionally be further smoothed using the function \code{\link{smoothDisplacementField}}. The smoothing kernels are  "Gauss","Laplace" and "Exponential". The displacement at point \code{x} will be the weighted displacment vectors of the k-closest displacement vectors. Be \code{d} the distance to a neightbouring point, the weight will be calculated as:
#' 
#' Gaussian:  \eqn{w(d) = exp(\frac{-d^2}{2\sigma^2})}{w(d) = exp(-d^2/2*sigma^2)}
#' 
#' Laplacian: \eqn{w(d) = exp(\frac{-d}{\sigma})}{w(d) = exp(-d/sigma)}
#'
#' Exponential: \eqn{w(d) = exp(\frac{-d}{2\sigma^2})}{w(d) = exp(-d/2*sigma^2)}
#' @seealso \code{\link{plot.DisplacementField}, \link{applyDisplacementField}, \link{smoothDisplacementField}}
#' @examples
#' require(Rvcg);require(Morpho)
#' data(dummyhead)
#' humoff <- meshOffset(dummyhead.mesh,offset=5)
#' dispfield <- createDisplacementField(dummyhead.mesh,humoff)
#' \dontrun{
#' ## this only runs with latest Rvcg build from master
#' highres <- vcgSubdivide(dummyhead.mesh)
#' ifield <- interpolateDisplacementField(dispfield,highres,threads=2,sigma = 50,k=500)
#' }
#' @export
interpolateDisplacementField <- function(dispfield, points, k=10, sigma=20,gamma=1,type=c("Gauss","Laplace","Exponential"), threads=parallel::detectCores()) {
    typeargs <- c("gauss","laplace","exponential")
    type <- match.arg(tolower(type[1]),typeargs)
    type <- match(type,typeargs)-1
    if (!inherits(dispfield,"DisplacementField"))
        stop("please enter displacementfield")
    if (inherits(points,"mesh3d"))
        points <- vert2points(points)
    if (is.vector(points))
        points <- t(points)
    ## get closest points on domain
    clost <- vcgKDtree(dispfield$domain,points,k=k)
    clost$index <- clost$index-1L
    tmp = .Call("smoothField",points, dispfield$DisplacementField,sigma,gamma,clost$index, clost$distance,threads,type)
    out <- list(domain=points,DisplacementField=tmp)
    class(out) <- "DisplacementField"
    ## if (smoothresult)
    ##    out <- smoothDisplacementField(out,sigma=sigma,k=k,threads = threads)
    return(out)
    
    
}

#' smooth a displacement field using Gaussian smoothing
#'
#' smooth a displacement field using Gaussian smoothing
#' @param dispfield displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @param k integer: number of k closest points to evaluate.
#' @param sigma sigma controls the weight of the neighbourhood by defining the standard-deviation for the gaussian smoothing
#' @param type kernel function for smoothing are "Gauss","Laplace" and "Exponential"
#' @param threads integer: number of threads to use for computing the interpolation.
#' @seealso \code{\link{interpolateDisplacementField}, \link{applyDisplacementField}, \link{plot.DisplacementField}}
#' @export
smoothDisplacementField <- function(dispfield,k=10,sigma=20,type=c("Gauss","Laplace","Exponential"),threads=1) {
    validDisplaceField(dispfield)
    typeargs <- c("gauss","laplace","exponential")
    type <- match.arg(tolower(type[1]),typeargs)
    type <- match(type,typeargs)-1
    gamma <- 1
    clost <- vcgKDtree(dispfield$domain,dispfield$domain,k=k)
    clost$index <- clost$index-1L
    tmp = .Call("smoothField",dispfield$domain, dispfield$DisplacementField,sigma,gamma,clost$index, clost$distance,threads,type)
    dispfield$DisplacementField <- tmp
    return(dispfield)
}



#' apply a discrete displacement field to a set of points/mesh in its domain
#'
#' apply a discrete displacement field to a set of points/mesh in its domain by applying the gaussian smoothed interpolation based of k closest neighbours
#' @param dispfield displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}} or 
#' @param points matrix or mesh3d at which to evaluate the interpolated displacement field
#' @param sigma sigma controls the weight of the neighbourhood by defining the standard-deviation for the gaussian smoothing
#' @param type kernel function for smoothing are "Gauss","Laplace" and "Exponential"
#' @param gamma dampening factor (displacement vectors will be divided by \code{gamma}
#' @param k integer: number of k closest points to evaluate.
#' @param threads integer: number of threads to use for computing the interpolation.
#' @note if points is identical to the domain of the displacement field, no interpolation will be performed.
#' @return returns the displaced version of points
#' @export
applyDisplacementField <- function(dispfield,points,k=10,sigma=20,type=c("Gauss","Laplace","Exponential"), gamma=1, threads=1) {
    validDisplaceField(dispfield)
    ## check if we need to interpolate at all
    if (!checkDispFieldDomain(dispfield,points))
        displacement <- interpolateDisplacementField(dispfield,points=points,k=k,sigma=sigma,type=type,gamma=gamma,threads=threads)
    else
        displacement <- dispfield
    updatePos <- displacement$domain+displacement$DisplacementField
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
#' @param x displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @param lwd width of the displacement vectors
#' @param color logical: if TRUE, displacement vectors will be colored according to a heatmap.
#' @param plot if FALSE only an object of class "DisplacementPlot" is returned without graphical output
#' @param ... additional arguments. Currently not used.
#' @return invisible object of class "DisplacementPlot" that can be rendered using \code{plot}
#' @seealso \code{\link{interpolateDisplacementField}, \link{applyDisplacementField}, \link{smoothDisplacementField}}
#' @importFrom Morpho plotNormals meshDist
#' @importFrom colorRamps blue2green2red
#' @importFrom rgl wire3d
#' @method plot DisplacementField
#' @export
plot.DisplacementField <- function(x,lwd=1,color=TRUE,plot=TRUE,...) {
     validDisplaceField(x)
     start <- x$domain
         end <- x$domain+x$DisplacementField
         vl <- nrow(start)
         dismesh <- list();class(dismesh) <- list("mesh3d","DisplacementPlot")
         dismesh$vb <- rbind(t(rbind(start,end)),1)
         dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
     if (!color) {
         if (plot)
          plot(dismesh,lit=FALSE,lwd=lwd)
     
     } else {
         dists <- apply(x$DisplacementField,1,norm,"2")
         ramp <- colorRamps::blue2green2red(19)
         from <- 0
         to <- ceiling(max(dists))
         colseq <- seq(from=from,to=to,length.out=20)
         coldif <- colseq[2]-colseq[1]
         distqual <- ceiling((dists/coldif)+1e-14)
         colorall <- ramp[distqual]
         dismesh$material$color <- rbind(colorall,colorall,colorall)
         dismesh$material$lit=FALSE
         if (plot)
             plot(dismesh,lwd=lwd)
     }
     invisible(dismesh)
     
}


#' plot the output of plot.DisplacementField
#'
#' plot the output of plot.DisplacementField
#' @param x displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @param lwd width of the displacement vectors
#' @param color logical: if TRUE, displacement vectors will be colored according to a heatmap.
#' @param ... additional arguments. Currently not used.
#' @seealso \code{\link{interpolateDisplacementField}, \link{applyDisplacementField}, \link{smoothDisplacementField}}
#' @examples
#' require(Rvcg)
#' data(dummyhead)
#' humoff <- meshOffset(dummyhead.mesh,offset=5)
#' dispfield <- createDisplacementField(dummyhead.mesh,humoff)
#' displot <- plot(dispfield,plot=FALSE)
#' plot(displot)
#' @method plot DisplacementPlot
#' @export
plot.DisplacementPlot <- function(x,lwd=1,color=TRUE,...) {
    nvert <- ncol(x$vb)
   
    if (!is.null(x$material$color)) {
        ptcol <- data.frame(as.vector(x$it),as.vector(x$material$color))
        ptcol <- unique(ptcol)
        ptcol <- ptcol[order(ptcol[,1]),2][1:(nvert/2)]
    } else {
        ptcol = 1
    }
    
    points3d(vert2points(x)[1:(nvert/2),],col=ptcol,size=lwd+4)
    wire3d(x,lwd=lwd)
    
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

## validates a displacement field
validDisplaceField <- function(x) {
    if (!inherits(x,"DisplacementField"))
        stop("not a valid displacement field")
    dim1 <- dim(x$domain)
    dim2 <- dim(x$DisplacementField)
    if (!isTRUE(all.equal(dim1,dim2)))
        stop("not a valid displacement field1")
}




