\name{colorTrans}
\alias{colorTrans}
\alias{mixColorMesh}

\title{
functions handling transfer of vertex colors of meshes.
}
\description{functions handling transfer of vertex colors of meshes.
}
\usage{
colorTrans(mesh1, mesh2, tol = 1)
mixColorMesh(mesh1, mesh2, alpha = 0.5)

}
\arguments{
  \item{mesh1}{triangular mesh of class "mesh3d".
%%     ~~Describe \code{mesh1} here~~
}
  \item{mesh2}{triangular mesh of class "mesh3d".
%%     ~~Describe \code{mesh2} here~~
}
\item{tol}{maximal distance of closest point included in color transfer.
}
\item{alpha}{numeric: 0 <= alpha <=1 weight of color associated with mesh1.

}
}
\details{colorTrans transfers the vertex color of mesh2 onto mesh1 by
  interpolating vertex colors of closest face.
  mixColorMesh mixes the vertex colors of two registered meshes
  (identical faces and corresponding vertices).
%%  ~~ If necessary, more details than the description above ~~
}
\value{colorTrans returns the mesh1 with updated vertex colors.

}

\author{ Stefan Schlager

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
