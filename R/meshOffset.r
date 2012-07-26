meshOffset <- function(mesh,offset)
{
  if (is.null(mesh$normals))
    mesh <- adnormals(mesh)

  mesh$vb[1:3,] <- mesh$vb[1:3,]+offset*mesh$normals[1:3,]
  invisible(mesh)
}
