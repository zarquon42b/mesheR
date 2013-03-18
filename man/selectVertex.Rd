\name{selectVertex}
\alias{selectVertex}
\alias{cutMesh}

\title{select and crop triangular surface meshes
}
\description{
select a subset of a surface mesh interactively.}
\usage{
selectVertex(mesh, col = 3, visible = TRUE, add = FALSE, render =
c("shade", "wire"), ...)
cutMesh(mesh, visible = TRUE, keep.selected = TRUE, col = 3, 
    add = FALSE, render = c("shade", "wire"), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"
%%     ~~Describe \code{mesh} here~~
}
  \item{col}{color to render the mesh.
%%     ~~Describe \code{col} here~~
}
  \item{visible}{select only vertices visible (from the present point of view)
%%     ~~Describe \code{visible} here~~
}
  \item{add}{logical: add the surface to an existing window.
%%     ~~Describe \code{add} here~~
}
  \item{render}{character: how to render the surface. Possible values
    are "shade" or "wire".
%%     ~~Describe \code{render} here~~
}
  \item{keep.selected}{logical: determines if the selected vertices or
    their complement are deleted.
%%     ~~Describe \code{\dots} here~~
} 
\item{\dots}{additional arguments passed to the rendering functions
    \code{\link{shade3d}} and \code{\link{wire3d}}
%%     ~~Describe \code{\dots} here~~
}
}
\details{select vertices or trim a triangular mesh by selecting a subset
  - either restricted to visible from the current point of view or by
  rectangular region. Visibility is determined by the fact that the
  vector from the viewpoint to the selected vertex does not intersect
  the mesh.
}
\value{
  selectVertex returns the indices of the selected vertices.

  cutMesh returns the trimmed mesh.

}

\author{
Stefan Schlager}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vcgPlyRead}}
}
\examples{
##Not run
data(nose)
## scale mesh to reasonable size
shortnose.mesh <- scalemesh(shortnose.mesh,500)
selection <- selectVertex(shortnose.mesh)
## End(Not run)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
