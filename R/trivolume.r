trivolume <- function(mesh1,mesh2)
  {
    #mesh2 <- conv2backf(mesh2)
    vb1 <- mesh1$vb
    vb2 <- mesh2$vb
    it <- mesh1$it; storage.mode(it) <- "integer"
    V <- 0; storage.mode(V) <- "double"
    out <- .Fortran("trianvol",vb1,vb2,it,ncol(vb1),ncol(it),V)
    return(out)
  }
