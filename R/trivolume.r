#' Calculate volume between two states of a mesh
#'
#' Calculate volume between two states of a mesh by dividing it into prisms
#'
#' @param mesh1 triangular mesh in first state
#' @param mesh2 triangular mesh in second state
#' @return Volume
#' @export trivolume
trivolume <- function(mesh1,mesh2)
  {
      if (missing(mesh2)) {
          mesh2 <- mesh1
      }
    #mesh2 <- conv2backf(mesh2)
    vb1 <- mesh1$vb
    vb2 <- mesh2$vb
    it <- mesh1$it-1; storage.mode(it) <- "integer"
    V <- 0; storage.mode(V) <- "double"
    out <- .Call("trianvol",vb1,vb2,it)
    return(out)
  }
