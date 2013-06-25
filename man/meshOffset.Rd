\name{meshOffset}
\alias{meshOffset}
\title{
inflate a mesh along its normals}
\description{
translate vertices of a triangular mesh along its normals}
\usage{
meshOffset(mesh, offset)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"
}
  \item{offset}{
distance to translate the vertices}
}

\value{returns modified mesh.
}

\author{
Stefan Schlager}



\examples{
data(nose)
offset <- meshOffset(shortnose.mesh,3)
shade3d(shortnose.mesh,col=3)
wire3d(offset)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
