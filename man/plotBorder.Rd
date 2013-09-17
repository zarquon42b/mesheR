\name{plotBorder}
\alias{plotBorder}
\title{plot all border edges of a triangular mesh}
\usage{
  plotBorder(mesh, col = 2, lwd = 2)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d".}

  \item{col}{color of lines.}

  \item{lwd}{line width}
}
\description{
  plot all border edges of a triangular mesh
}
\examples{
require(rgl)
require(Morpho)
data(nose)
plotBorder(shortnose.mesh)
wire3d(shortnose.mesh,col=3)
}
\author{
  Stefan Schlager
}
\keyword{~kwd1}
\keyword{~kwd2}

