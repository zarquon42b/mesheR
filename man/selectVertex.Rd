\name{selectVertex}
\alias{cutMesh}
\alias{selectVertex}
\title{select and crop triangular surface meshes}
\usage{
  selectVertex(mesh, col = 3, visible = TRUE, add = FALSE,
    render = c("shade", "wire"), ...)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"}

  \item{col}{color to render the mesh.}

  \item{visible}{select only vertices visible (from the
  present point of view)}

  \item{add}{logical: add the surface to an existing
  window.}

  \item{render}{character: how to render the surface.
  Possible values are "shade" or "wire".}

  \item{keep.selected}{logical: determines if the selected
  vertices or their complement are deleted.}

  \item{\dots}{additional arguments passed to the rendering
  functions \code{\link{shade3d}} and \code{\link{wire3d}}}
}
\value{
  selectVertex returns the indices of the selected
  vertices.

  cutMesh returns the trimmed mesh.
}
\description{
  select a subset of a surface mesh interactively.
}
\details{
  select vertices or trim a triangular mesh by selecting a
  subset - either restricted to visible from the current
  point of view or by rectangular region. Visibility is
  determined by the fact that the vector from the viewpoint
  to the selected vertex does not intersect the mesh.
}
\examples{
\dontrun{
data(nose)
selection <- selectVertex(shortnose.mesh)
}
}
\author{
  Stefan Schlager
}
\seealso{
  \code{\link{vcgPlyRead}}
}
\keyword{~kwd1}
\keyword{~kwd2}

