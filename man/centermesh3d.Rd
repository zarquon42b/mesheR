\name{centermesh3d}
\alias{centermesh3d}
\title{move the centroid of a mesh to a given coordinate}
\usage{
  centermesh3d(x, center = NULL)
}
\arguments{
  \item{x}{triangular mesh of class 'mesh3d'}

  \item{center}{coordinate where to translate the mesh, if
  center=NULL the centroid will be moved to the origin}
}
\value{
  returns the translated mesh
}
\description{
  move the centroid of a mesh to a given coordinate
}
\examples{
data(humface)
wire3d(humface, col=3)
humcenter <- centermesh3d(humface)
#view translated mesh
wire3d(humcenter, col=2)
}
