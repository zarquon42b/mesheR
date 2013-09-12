#' make a warpmovie with multiple registered meshes
#'
#' make a warpmovie with multiple registered meshes - this is a pimped version of warpmovie3d from package 'Morpho'
#'
#' @param \dots registered meshes
#' @param n amount of intermediate images between two meshes
#' @param col mesh color
#' @param folder character: output folder for created images (optional)
#' @param movie character: name of the output files
#' @param add logical: if TRUE, the movie will be added to the focussed rgl-windows.
#' @param close logical: if TRUE, the rgl window will be closed when finished. width and 200 the height of the image.
#' @param whichcolor select which mesh's vertexcolors to use. 
#' 
#' @param countbegin integer: number to start image sequence. 
#' @param ask logical: if TRUE, the viewpoint can be selected manually. 
#' @importFrom rgl points3d open3d rgl.pop rgl.snapshot rgl.close shade3d
#' @export warpmovieMulti
warpmovieMulti <- function(..., n, col="green", folder=NULL, movie="warpmovie",add=FALSE, close=TRUE, countbegin=0, ask=TRUE, whichcolor=1)
{	
    args <- list(...)
    argc <- length(args)

    if(!is.null(folder)) {
        if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") {
            folder<-paste(folder,"/",sep="")
            dir.create(folder,showWarnings=F)
            movie <- paste(folder,movie,sep="")
        }
    }
    if (!add)
        open3d()
    
    ## get bbox
    ##bbox <- apply(rbind(vert2points(x),vert2points(y)),2,range)
    bbox0 <- lapply(args,function(x) x <- apply(vert2points(x),2,range))
    bbox <- vert2points(args[[1]])
    for (i in 2:length(args))
        bbox <- rbind(bbox,vert2points(args[[i]]))
    bbox <- apply(bbox,2,range)
    bbox <- expand.grid(bbox[,1],bbox[,2],bbox[,3])
    points3d(bbox,col="white",alpha=0)

    for (j in 1:(argc-1)) {
        x <- args[[j]]
        y <- args[[j+1]]
        if (!is.null(args[[whichcolor]]$material$color)) {
            y$material$color <- args[[whichcolor]]$material$color
            x$material$color <- args[[whichcolor]]$material$color
        }
        
        for (i in 0:n) {
            mesh<-x
            mesh$vb[1:3,]<-(i/n)*y$vb[1:3,]+(1-(i/n))*x$vb[1:3,]
            mesh <- adnormals(mesh)
            a <- shade3d(mesh,col=col,specular=1)
            if (i ==0 && j == 1 && ask==TRUE)
                readline("please select view and press return\n")
                                     
            filename <- sprintf("%s%04d.png", movie, countbegin+i)
            rgl.snapshot(filename,fmt="png")
            rgl.pop("shapes",id=a)
        }
        countbegin <- countbegin+n
    }
    
    if (FALSE) {#palindrome) ## go the other way ##
        for (i in 1:(n-1)) {
            mesh<-x
            mesh$vb[1:3,]<-(i/n)*x$vb[1:3,]+(1-(i/n))*y$vb[1:3,]
            mesh <- adnormals(mesh)
            a <- shade3d(mesh,col=col, shadeopts)
            filename <- sprintf("%s%04d.png", movie, countbegin+i+n)
            rgl.snapshot(filename,fmt="png")
            rgl.pop("shapes",id=a)	
        }
    }
    if (close)
        rgl.close()
}
