\name{rmInternals}
\alias{rmInternals}
\title{remove vertices within a convex surface mesh
}
\description{
remove vertices within a convex surface mesh by checking if they
intersect the mesh along their normals.
}
\usage{
rmInternals(mesh, mindist = NULL,explode=FALSE,center=NULL)
}

\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"
%%     ~~Describe \code{mesh} here~~
}
  \item{mindist}{numeric: restrict to those vertices whose intersection point is
  further away as mindist.
%%     ~~Describe \code{mindist} here~~
}  
\item{explode}{logical: instead of surface normals use the direction
  outward from a center of the mesh. Per default, the center of the
  meshes bounding box is used.
%%     ~~Describe \code{mindist} here~~
}
\item{center}{numeric vector of length 3, used in the case
  "explode=TRUE", if another center is required.
%%     ~~Describe \code{mindist} here~~
}
}
\details{This only works for convex meshes, for concave substructures
  this might result in unwanted removals.
}
\value{returns a cleaned mesh

}

\author{Stefan Schlager
}

\seealso{
\code{\link{rmVertex}}
}
%\examples{

%  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
