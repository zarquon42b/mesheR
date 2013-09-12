\name{glVisible}
\alias{glVisible}
\title{determine visibility of mesh's vertices}
\usage{
  glVisible(mesh, offset = 0.001)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d". Must be
  currently rendered in an rgl window.}

  \item{offset}{initial offset to move vertex slightly away
  from the surface.}
}
\value{
  returns logical vector, assigning TRUE/FALSE to each
  vertex of a mesh.
}
\description{
  Determine which vertices of a rendered mesh are visible
  from the current viewpoint
}
\examples{
require(rgl)
data(nose)
shade3d(shortnose.mesh,col=3)
visi <- glVisible(shortnose.mesh)
points3d(vert2points(shortnose.mesh)[which(visi),])
}
\author{
  Stefan Schlager
}
\seealso{
  \code{\link{selectVertex}}, \code{\link{cutMesh}}
}
\keyword{~kwd1}
\keyword{~kwd2}

