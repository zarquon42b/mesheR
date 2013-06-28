\name{colorTrans}
\alias{colorTrans}
\alias{mixColorMesh}
\title{functions handling transfer of vertex colors of meshes.}
\usage{
  colorTrans(mesh1, mesh2, tol = 1)
}
\arguments{
  \item{mesh1}{triangular mesh of class "mesh3d".}

  \item{mesh2}{triangular mesh of class "mesh3d".}

  \item{tol}{maximal distance of closest point included in
  color transfer.}

  \item{alpha}{numeric: 0 <= alpha <=1 weight of color
  associated with mesh1.}
}
\value{
  colorTrans returns the mesh1 with updated vertex colors.
}
\description{
  colorTrans transfers the vertex color of mesh2 onto mesh1
  by interpolating vertex colors of closest face.
  mixColorMesh mixes the vertex colors of two registered
  meshes (identical faces and corresponding vertices).
}
\examples{
data(nose)
redmesh <- shortnose.mesh
redmesh$material$color <- matrix("#FF0000",dim(shortnose.mesh$it))
bluemesh <- shortnose.mesh
bluemesh$material$color <- matrix("#0000FF",dim(shortnose.mesh$it))
mixmesh <- mixColorMesh(bluemesh,redmesh)
shade3d(mixmesh)
}
\author{
  Stefan Schlager
}
\keyword{~kwd1}
\keyword{~kwd2}

