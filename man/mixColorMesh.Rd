\name{mixColorMesh}
\alias{mixColorMesh}
\title{mix vertex colors of two registered meshes.}
\usage{
  mixColorMesh(mesh1, mesh2, alpha = 0.5)
}
\arguments{
  \item{mesh1}{triangular mesh of class "mesh3d".}

  \item{mesh2}{triangular mesh of class "mesh3d".}

  \item{alpha}{numeric: 0 <= alpha <=1 weight of color
  associated with mesh1.}
}
\value{
  returns mesh1 with mixed vertex colors.
}
\description{
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

