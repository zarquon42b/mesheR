\name{meshOffset}
\alias{meshOffset}
\title{
inflate a mesh}
\description{
inflate or deflate a mesh along its normals}
\usage{
meshOffset(mesh, offset)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"
}
  \item{offset}{numeric: amount of in/deflation (positive values go
    outward, while negative values go inward)
}
}

\value{
object of class "mesh3d"
}

\author{
Stefan Schlager}

}
\examples{
data(nose)
noseinf <- meshOffset(shortnose.mesh,0.001)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
