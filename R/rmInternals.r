rmInternals <- function(mesh,mindist=NULL)
  {
    meshoff <- meshOffset(mesh,0.05)
    check <- vcgIntersect(meshoff,mesh)
    if (is.null(mindist))
      intern <- which(as.logical(check$quality))
    else
      intern <- which((check$quality) * (check$distance > mindist))

    out <- rmVertex(mesh,intern)
    cat(paste("removed",length(intern),"vertices\n"))
    invisible(out)
  }

    
