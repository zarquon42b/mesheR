glView <- function()
  {
    out <- c(0,0,0)
    info <- rgl.projection()
    mdl <- info$model
    out[1] <- -(mdl[1] * mdl[13] + mdl[2] * mdl[14] + mdl[3] * mdl[15])
    out[2] <-  -(mdl[5] * mdl[13] + mdl[6] * mdl[14] + mdl[7] * mdl[15])
    out[3] <- -(mdl[9] * mdl[13] + mdl[10] * mdl[14] + mdl[11] * mdl[15])
    return(out)
  }
glVisible <- function(mesh)
{
  mesh <- adnormals(mesh)
  mesh0 <- meshOffset(mesh,1e-5)
  viewpoint <- c(glView(),0)
  normals <- viewpoint-mesh0$vb
  mesh0$normals <- normals
  tmp <- as.logical(vcgIntersect(mesh0,mesh)$quality)
  out <- tmp
  out[tmp] <- FALSE
  out[!tmp] <- TRUE
  return(out)
}

    
