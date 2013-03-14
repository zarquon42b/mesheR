rmInternals <- function(mesh,mindist=NULL,explode=FALSE,center=NULL)
  {

    
    if (explode)
      {
        if (is.null(center))
          center <- apply(meshcube(mesh),2,mean)
        mesh$normals[1:3,] <- mesh$vb[1:3,]-center
        mesh$normals[1:3,] <- apply(mesh$normals[1:3,],2,function(x){x <- x/(sqrt(sum(x^2)))})
      }
    meshoff <- meshOffset(mesh,0.005)
    check <- vcgIntersect(meshoff,mesh)
    if (is.null(mindist))
      intern <- which(as.logical(check$quality))
    else
      intern <- which(((check$quality) * (check$distance > mindist))==1)

    out <- rmVertex(mesh,intern)
    cat(paste("removed",length(intern),"vertices\n"))
    invisible(out)
  }

    
