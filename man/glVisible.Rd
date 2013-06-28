\name{glVisible}
\alias{glVisible}
\title{determine visibility of mesh's vertices}
\usage{
  glVisible(mesh)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d". Must be
  currently rendered in an rgl window.}
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
data(nose)
shade3d(shortnose.mesh,col=3)
visi <- glVisible(shortnose.mesh)
points3d(vert2points(shortnose.mesh)[which(visi),])
}
\author{
  Stefan Schlager
}
\keyword{~kwd1}
\keyword{~kwd2}

