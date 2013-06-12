\name{plotBorder}
\alias{plotBorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  plot all border edges of a triangular mesh}
\description{plot all border edges of a triangular mesh
}
\usage{
plotEdges(mesh, col = 2, lwd = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d".
    %%     ~~Describe \code{mesh} here~~
  }
  \item{col}{color of lines.
    %%     ~~Describe \code{col} here~~
  }
  \item{lwd}{line width
    %%     ~~Describe \code{lwd} here~~
  }
}


\author{
  Stefan Schlager}

\examples{
data(nose)
plotBorder(shortnose.mesh)
wire3d(shortnose.mesh,col=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
