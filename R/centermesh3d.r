#' @export centermesh3d
centermesh3d <- function(x,center=NULL)
  {

    bary <- apply(vert2points(x),2,mean)
    if (is.null(center))
     center <- c(0,0,0)
     
       
    x$vb[1:3,] <- x$vb[1:3,]- (bary-center)
    invisible(x)
  }
 
