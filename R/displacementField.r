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
#' @param sigma kernel bandwidth used for smoothing. For all kernels except B-spline, sigma controls the importance of the neighbourhood by defining the bandwidth of the smoothing kernel. For B-spline it defines the support (the higher, the "wobblier" the deformation field can become.
#' @param gamma dampening factor (displacement vectors will be divided by \code{gamma}
#' @param type kernel function for smoothing are "Gauss","Laplace", "Exponential", "Bspline" and "TPS" (or any abbreviation thereof).
#' @param subsample integer: amount to subsample the field in case of type="TPS"
#' @param lambda smoothing factor for TPS
#' @param threads integer: number of threads to use for computing the interpolation.
#' @param ... additional parameters - currentyl unused.
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
#' ifield <- interpolateDisplacementField(dispfield,highres,threads=2,sigma = 10,k=50)
#' }
#' @importFrom Morpho tps3d
#' @export
interpolateDisplacementField <- function(dispfield, points, k=10, sigma=20,gamma=1,type=c("Gauss","Laplace","Exponential","Bspline","TPS"),subsample=2000, lambda=1e-8,threads=0,...) {
    validDisplaceField(dispfield)
    typeargs <- c("gauss","laplace","exponential","bspline","tps")
    type <- match.arg(tolower(type[1]),typeargs)
    type <- match(type,typeargs)-1
    if (!inherits(dispfield,"DisplacementField"))
        stop("please enter displacementfield")
    if (inherits(points,"mesh3d"))
        points <- vert2points(points)
    if (is.vector(points))
        points <- t(points)
    ## get closest points on domain
    clost <- vcgKDtree(dispfield$domain,points,k=k,threads=threads)
    clost$index <- clost$index-1L
    if (type %in% 0:3) {
        tmp = .Call("smoothField",points, dispfield$domain, dispfield$DisplacementField,sigma,gamma,clost$index, clost$distance,iterations=1,threads,type)
    } else {
        subind <- fastKmeans(dispfield$domain,k=subsample,threads=threads)
        tar <- dispfield$domain[subind$selected,,drop=FALSE]+dispfield$DisplacementField[subind$selected,,drop=FALSE]
        tmp <- tps3d(points,dispfield$domain[subind$selected,,drop=FALSE],tar,threads = threads,lambda=lambda)
        tmp <- tmp-points
    }
    
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
#' @param sigma kernel bandwidth used for smoothing. For all kernels except B-spline, sigma controls the importance of the neighbourhood by defining the bandwidth of the smoothing kernel. For B-spline it defines the support (the higher, the "wobblier" the deformation field can become.
#' @param type kernel function for smoothing are "Gauss","Laplace", "Exponential", "Bspline" and "TPS" (or any abbreviation thereof).
#' @param iterations number of iterations to run
#' @param subsample integer: amount to subsample the field in case of type="TPS"
#' @param lambda smoothing factor for TPS
#' @param threads integer: number of threads to use for computing the interpolation.
#' @seealso \code{\link{interpolateDisplacementField}, \link{applyDisplacementField}, \link{plot.DisplacementField}}
#' @export
smoothDisplacementField <- function(dispfield,k=10,sigma=20,type=c("Gauss","Laplace","Exponential","Bspline","TPS"),iterations=1,subsample=2000, lambda=1e-8,threads=0) {
    validDisplaceField(dispfield)
    typeargs <- c("gauss","laplace","exponential","bspline","tps")
    type <- match.arg(tolower(type[1]),typeargs)
    type <- match(type,typeargs)-1
    gamma <- 1
    clost <- vcgKDtree(dispfield$domain,dispfield$domain,k=k,threads=threads)
    clost$index <- clost$index-1L
    if (type %in% 0:3) {
        tmp = .Call("smoothField",dispfield$domain, dispfield$domain, dispfield$DisplacementField,sigma,gamma,clost$index, clost$distance,iterations,threads,type)
    } else {
        subind <- Morpho::fastKmeans(dispfield$domain,k=subsample,threads=threads)
        tar <- dispfield$domain[subind$selected,,drop=FALSE]+dispfield$DisplacementField[subind$selected,,drop=FALSE]
        tmp <- tps3d(dispfield$domain,dispfield$domain[subind$selected,,drop=FALSE],tar,threads = threads,lambda=lambda)
        tmp <- tmp-dispfield$domain
    }
    dispfield$DisplacementField <- tmp
    return(dispfield)
}



#' apply a discrete displacement field to a set of points/mesh in its domain
#'
#' apply a discrete displacement field to a set of points/mesh in its domain by applying the gaussian smoothed interpolation based of k closest neighbours
#' @param dispfield displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}} or 
#' @param points matrix or mesh3d at which to evaluate the interpolated displacement field
#' @param sigma kernel bandwidth used for smoothing. For all kernels except B-spline, sigma controls the importance of the neighbourhood by defining the bandwidth of the smoothing kernel. For B-spline it defines the support (the higher, the "wobblier" the deformation field can become.
#' @param type kernel function for smoothing are "Gauss","Laplace", "Exponential" and "Bspline" (or any abbreviation thereof).
#' @param gamma dampening factor (displacement vectors will be divided by \code{gamma}
#' @param k integer: number of k closest points to evaluate.
#' @param lambda smoothing factor for TP
#' @param threads integer: number of threads to use for computing the interpolation.
#' @note if points is identical to the domain of the displacement field, no interpolation will be performed.
#' @return returns the displaced version of points
#' @export
applyDisplacementField <- function(dispfield,points,k=10,sigma=20,type=c("Gauss","Laplace","Exponential","TPS"), gamma=1,lambda=1e-8, threads=1) {
    validDisplaceField(dispfield)
    ## check if we need to interpolate at all
    if (!checkDispFieldDomain(dispfield,points,threads=threads)) {
        message("dispfield domain and points not identical: interpolating...")
        displacement <- interpolateDisplacementField(dispfield,points=points,k=k,sigma=sigma,type=type,gamma=gamma,lambda=lambda,threads=threads)
    } else
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
plot.DisplacementField <- function(x,lwd=1,color=TRUE,plot=TRUE,size=NULL,...) {
    validDisplaceField(x)
    start <- x$domain
    end <- x$domain+x$DisplacementField
    vl <- nrow(start)
    dismesh <- list();class(dismesh) <- list("mesh3d","DisplacementPlot")
    dismesh$vb <- rbind(t(rbind(start,end)),1)
    dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
    if (!color) {
        if (plot)
            plot(dismesh,lit=FALSE,lwd=lwd,size=size,...)
        
    } else {
        dists <- sqrt(rowSums(x$DisplacementField^2))
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
            plot(dismesh,lwd=lwd,size=size,...)
    }
    invisible(dismesh)
    
}


#' plot the output of plot.DisplacementField
#'
#' plot the output of plot.DisplacementField
#' @param x displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @param lwd width of the displacement vectors
#' @param color logical: if TRUE, displacement vectors will be colored according to a heatmap.
#' @param size size of grid points
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
plot.DisplacementPlot <- function(x,lwd=1,color=TRUE,size=NULL,...) {
    nvert <- ncol(x$vb)
    if (!is.null(x$material$color)) {
        ptcol <- data.frame(as.vector(x$it),as.vector(x$material$color))
        ptcol <- unique(ptcol)
        ptcol <- ptcol[order(ptcol[,1]),2][1:(nvert/2)]
    } else {
        ptcol = 1
    }
    if (is.null(size))
        size <- lwd+4
    points3d(vert2points(x)[1:(nvert/2),],col=ptcol,size=size)
    wire3d(x,lwd=lwd)
    
}


## checks whether a set of points is identical to the domain of a displacement field
checkDispFieldDomain <- function(dispfield,points,tol=1e-12,threads=1) {
    if (inherits(points,"mesh3d"))
        points <- vert2points(points)
    ndisp <- nrow(dispfield$domain)
    npoints <- nrow(points)
    if (ndisp != npoints)
        return(FALSE)
    else {
        chk <- vcgKDtree(dispfield$domain,points,k=1,threads=threads)
        if (!isTRUE(all.equal(as.vector(chk$index),seq.int(ndisp),check.attributes = FALSE)))
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



#' convert a discrete irregular displacement field to a displacement grid
#'
#' convert a discrete irregular displacement field to a displacement grid
#' @param x field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}}
#' @param spacing spacing between grid nodes
#' @param margin percentage of offset around \code{x}
#' @param IJK2RAS image to coordinate transform
#' @param invert if TRUE, the displacement field will be inverted
#' @param ... additional arguments passed to \code{\link{interpolateDisplacementField}} to choose interpolation etc..
#' @return an object of class DisplacementField with additional attributes
#' @export
displacementField2Grid <- function(x,spacing=rep(2,3),margin=0.05,IJK2RAS=diag(c(-1,-1,1,1)),invert=FALSE,...) {
    if (!inherits(x,"DisplacementField"))
        stop("x must be a displacement field")
    if (invert)
        x <- invertDisplacementField(x)
    xtrans <- lapply(x,applyTransform,IJK2RAS)
    class(xtrans) <- class(x)
    x <- xtrans
    pts <- (x$domain)
    ranges <- apply(pts,2,range)
    ranges <- apply(ranges,2,extendrange,f=margin)
    grid <- lapply(1:3,function(x) seq(from=ranges[1,x],to=ranges[2,x],by=spacing[x]))
    mygrid <- as.matrix(expand.grid(grid[[1]],grid[[2]],grid[[3]]))
    arrdims <- sapply(grid,length)
    myarr <- array(NA,dim=arrdims)
    indices <- as.matrix(expand.grid(1L:arrdims[1],1L:arrdims[2],1L:arrdims[3]))
    dfgrid <- interpolateDisplacementField(x,mygrid,...)
    
    gridattributes <- list(indices=indices,origin=mygrid[1,],arraydim=arrdims,spacing=spacing)
    attributes(dfgrid) <- append(attributes(dfgrid),gridattributes)
    class(dfgrid) <- c("DisplacementGrid","DisplacementField")
    return(dfgrid)
}


displacementGrid2Transform <- function(x) {
    if (!require(SimpleITK))
        stop("this function requires SimpleITK to be installed")
    if (!inherits(x,"DisplacementGrid"))
        stop("x must be of class DisplacementGrid")
    attribs <- attributes(x)
    img <- SimpleITK::Image(attribs$arraydim[1],attribs$arraydim[2],attribs$arraydim[3],"sitkVectorFloat64")
    inds <- attribs$indices-1L
    for(i in 1:(nrow(attribs$indices))) {
        tmp=as.numeric(x$DisplacementField[i,])
        img$SetPixelAsVectorFloat64(inds[i,],tmp)
    }
    img$SetSpacing(attribs$spacing)
    img$SetOrigin(attribs$origin)
    return(img)
}


