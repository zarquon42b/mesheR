#' make a warpmovie with multiple registered meshes
#'
#' make a warpmovie with multiple registered meshes - this is a pimped version of warpmovie3d from package 'Morpho'
#'
#' @param \dots registered meshes or a list containing registered meshes.
#' @param n amount of intermediate images between two meshes
#' @param col mesh color - overrides options \code{whichcolor} and \code{mixcolor}
#' @param folder character: output folder for created images (optional)
#' @param movie character: name of the output files
#' @param add logical: if TRUE, the movie will be added to the focussed rgl-windows.
#' @param close logical: if TRUE, the rgl window will be closed when finished. width and 200 the height of the image.
#' @param whichcolor select which mesh's vertexcolors to use.
#' @param align logical: align all meshes to the first one
#' @param scale if \code{align=TRUE} this controls if scaling is allowed.
#' 
#' @param countbegin integer: number to start image sequence. 
#' @param ask logical: if TRUE, the viewpoint can be selected manually.
#' @param mixcolor logical: if available existing colors are subsequently mixed using \code{\link{mixColorMesh}}.
#' @param shade select shade: "s"=shade3d, "w"= wire3d, "b"=both
#' @examples
#' require(Morpho)
#' data(nose)
#' redmesh <- shortnose.mesh
#' redmesh$material$color <- matrix("#FF0000",dim(shortnose.mesh$it))
#' bluemesh <- shortnose.mesh
#' bluemesh$material$color <- matrix("#0000FF",dim(shortnose.mesh$it))
#' \dontrun{
#' warpmovieMulti(bluemesh, redmesh, n=15)
#' }
#' @importFrom rgl open3d points3d shade3d rgl.snapshot rgl.pop rgl.close
#' @export
warpmovieMulti <- function(..., n, col=NULL, folder=NULL, movie="warpmovie",add=FALSE, close=TRUE, countbegin=0, ask=TRUE, whichcolor=NULL, align=TRUE, scale=FALSE, mixcolor=TRUE, shade=c("s","w","b")) UseMethod("warpmovieMulti")

#' @export
warpmovieMulti.default <- function(..., n, col=NULL, folder=NULL, movie="warpmovie",add=FALSE, close=TRUE, countbegin=0, ask=TRUE, whichcolor=NULL, align=TRUE, scale=FALSE, mixcolor=TRUE, shade=c("s","w","b")) {
    args <- list(...)
    warpmovieMulti(args,n=n, col=col, folder=folder, movie=movie,add=add, close=close, countbegin=countbegin, ask=ask, whichcolor=whichcolor, align=align, scale=scale, mixcolor=mixcolor, shade=shade)
}

#' @export    
warpmovieMulti.list <- function(..., n, col=NULL, folder=NULL, movie="warpmovie",add=FALSE, close=TRUE, countbegin=0, ask=TRUE, whichcolor=NULL, align=TRUE, scale=FALSE, mixcolor=TRUE, shade=c("s","w","b"))
{	
    args <- (...)
    argc <- length(args)
    if (argc < 2)
            stop("at least two arguments needed")

    npics <- nchar(as.character(argc*n))+1
    ndec <- paste0("%s%0",npics,"d.png")

    if(!is.null(folder)) {
        if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") 
            folder<-paste(folder,"/",sep="")
        dir.create(folder,showWarnings=F)
        movie <- paste(folder,movie,sep="")
        
    }
    if (!is.null(col)) {
        args[[1]] <- colorMesh(args[[1]],col)
        mixcolor <- FALSE
        whichcolor <- 1
    }  else
        col <- "white"
    
    if (!add)
        open3d()
    back <- front <- "filled"
    if (shade[1] == "w")
        back <- front <- "lines"
    ## get bbox
    ##bbox <- apply(rbind(vert2points(x),vert2points(y)),2,range)
    bbox0 <- lapply(args,function(x) x <- apply(vert2points(x),2,range))
    bbox <- vert2points(args[[1]])
    for (i in 2:length(args)) {
        if (align)
            args[[i]] <- rotmesh.onto(args[[i]], vert2points(args[[i]]), vert2points(args[[1]]),scale=scale, adnormals=TRUE)$mesh
        bbox <- rbind(bbox,vert2points(args[[i]]))
    }
    bbox <- apply(bbox,2,range)
    bbox <- expand.grid(bbox[,1],bbox[,2],bbox[,3])
    points3d(bbox,col="white",alpha=0)

    for (j in 1:(argc-1)) {
        x <- args[[j]]
        y <- args[[j+1]]
        if (!is.null(whichcolor)) {
            if (!is.null(args[[whichcolor]]$material$color)) {
                y$material$color <- args[[whichcolor]]$material$color
                x$material$color <- args[[whichcolor]]$material$color
            }
        }
                
        for (i in 0:n) {
            mesh <- x
            mesh$vb[1:3,]<-(i/n)*y$vb[1:3,]+(1-(i/n))*x$vb[1:3,]
            mesh <- vcgUpdateNormals(mesh)
            if (mixcolor && is.null(whichcolor)) {
                if (is.null(mesh$material$color)) 
                    mesh$material$color <- matrix("#FFFFFF",nrow(mesh$it),ncol(mesh$it))
                if (is.null(y$material$color)) 
                    y$material$color <- matrix("#FFFFFF",nrow(y$it),ncol(y$it))
                mesh <- mixColorMesh(mesh,y, alpha=((i/n)))
            }
                    
                
            a <- shade3d(mesh,col=col,specular=1,back=back,front=front)
            if (shade[1] == "b"){
                tmesh <- mesh
                tmesh$material$color <- NULL
                a <- append(a,wire3d(tmesh,col=1,lwd=2,specular=1))
            }
            if (i ==0 && j == 1 && ask==TRUE && interactive())
                readline("please select view and press return\n")
                                     
            filename <- sprintf(ndec, movie, countbegin+i)
            rgl.snapshot(filename,fmt="png")
            rgl.pop("shapes",id=a)
        }
        countbegin <- countbegin+n
    }
    
    if (FALSE) {#palindrome) ## go the other way ##
        for (i in 1:(n-1)) {
            mesh<-x
            mesh$vb[1:3,]<-(i/n)*x$vb[1:3,]+(1-(i/n))*y$vb[1:3,]
            mesh <- vcgUpdateNormals(mesh)
            a <- shade3d(mesh,col=col, shadeopts)
            filename <- sprintf("%s%04d.png", movie, countbegin+i+n)
            rgl.snapshot(filename,fmt="png")
            rgl.pop("shapes",id=a)	
        }
    }
    if (close)
        rgl.close()
}
