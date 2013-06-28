\name{AmbergDeformSpam}
\alias{AmbergDeformSpam}
\title{Deform triangular mesh based on correspondences}
\usage{
  AmbergDeformSpam(mesh, lm1, lm2, k0 = 1, lambda = 1,
    S = NULL)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d". No loose
  vertices, edges and degenerated faces are allowed.}

  \item{lm1}{m x 3 matrix containing correspondences on
  "mesh"}

  \item{lm2}{m x 3 matrix containing target
  correspondences}

  \item{k0}{integer: parameter regularizing face normal
  distortion.}

  \item{lambda}{numeric: parameter regularizing faces's
  distortion.}

  \item{S}{optional: object from function createS from
  previous calculation.}
}
\value{
  \item{mesh}{deformed mesh} \item{Jn}{Jacobi submatrix Jn}
  \item{Jc}{Jacobi submatrix Jc} \item{J }{Jacobian matrix}
  \item{H }{Hessian of J, class "spam"}
  \item{Hchol}{Cholesky decomposition of H; class"spam"}
}
\description{
  Perform smooth deformation of a triangular mesh,
  minimizing per-face distortions.
}
\details{
  Perform smooth deformation of a triangular mesh,
  minimizing per-face distortions.No loose vertices, edges
  and degenerated faces are allowed, as they lead to
  singular equation system.
}
\author{
  Stefan Schlager
}
\references{
  Amberg, B. 2011. Editing faces in videos, University of
  Basel.
}
\seealso{
  \code{\link{gaussMatch}}
}
\keyword{~kwd1}
\keyword{~kwd2}

