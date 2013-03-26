mesh2sparse <- function(mesh)
  {
    nvb <- dim(mesh$vb)[2]
    
    out <- sparseMatrix(j=1:(nvb*4),i=rbind(1:nvb,1:nvb,1:nvb,1:nvb),x=as.vector(mesh$vb[,]))
    invisible(out)
  }
