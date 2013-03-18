\name{glVisible}
\alias{glVisible}
\title{determine visibility of mesh's vertices
  
}
\description{
Determine which vertices of a rendered mesh are visible from the current
viewpoint}
\usage{
glVisible(mesh)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d". Must be currently
    rendered in an rgl window.
}
}


\value{returns logical vector, assigning TRUE/FALSE to each vertex of a mesh.
}


\author{
Stefan Schlager}

\examples{
data(nose)
##scale mesh to reasonable size
shortnose.scale <- scalemesh(shortnose.mesh,500)
shade3d(shortnose.scale,col=3)
visi <- glVisible(shortnose.scale)
points3d(vert2points(shortnose.scale)[which(visi),])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
