\name{colorTrans}
\alias{colorTrans}
\title{Transfers vertex colors between meshes.}
\usage{
  colorTrans(mesh1, mesh2, tol = 1, cores = 1)
}
\arguments{
  \item{mesh1}{triangular mesh of class "mesh3d".}

  \item{mesh2}{triangular mesh of class "mesh3d".}

  \item{tol}{maximal distance of closest point included in
  color transfer.}

  \item{cores}{maximum cores for \code{\link{closemeshKD}}
  to use.}
}
\value{
  returns the mesh1 with transfered vertex colors.
}
\description{
  colorTrans transfers the vertex color of mesh2 onto mesh1
  by interpolating vertex colors of closest face.
}
\examples{
data(nose)
redmesh <- shortnose.mesh
#color mesh red
redmesh$material$color <- matrix("#FF0000",dim(shortnose.mesh$it))
# mesh without colors
nocolmesh <- shortnose.mesh
mixmesh <- colorTrans(nocolmesh, redmesh)
shade3d(mixmesh)
}
\author{
  Stefan Schlager
}
\keyword{~kwd1}
\keyword{~kwd2}

