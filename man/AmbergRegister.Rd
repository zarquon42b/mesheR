\name{AmbergRegister}
\alias{AmbergRegister}
\title{Register two triangular meshes based on smooth deformation.}
\usage{
  AmbergRegister(mesh1, mesh2, lm1 = NULL, lm2 = NULL,
    k = 1, lambda = 1, iterations = 15, rho = pi/2,
    dist = 2, border = FALSE, smooth = TRUE, smoothit = 1,
    smoothtype = "t", tol = 1e-04, useiter = TRUE,
    minclost = 50, distinc = 1, scale = TRUE, icp = NULL)
}
\arguments{
  \item{mesh1}{reference mesh: triangular mesh of class
  "mesh3d". No loose vertices, edges and degenerated faces
  are allowed.}

  \item{mesh2}{target mesh: triangular mesh of class
  "mesh3d".}

  \item{lm1}{m x 3 matrix containing correspondences on
  "mesh1"}

  \item{lm2}{m x 3 matrix containing target correspondences
  on "mesh2"}

  \item{k}{integer: parameter regularizing face normal
  distortion. Can be vector of length(iterations) or single
  value.}

  \item{lambda}{numeric: parameter regularizing faces's
  distortion. Can be vector of length(iterations) or single
  value.}

  \item{iterations}{integer: number of iterations to run.}

  \item{rho}{numeric: 0 < rho < 2*pi tolerance of normal
  deviation between reference vertices and corresponding
  closest points on target suface.}

  \item{dist}{numeric: tolerance of maximal distance
  between reference vertices and corresponding closest
  points on target suface.}

  \item{border}{logical: if FALSE, hits on border faces are
  ignored (reduces distortion)}

  \item{smooth}{logical: if TRUE after each iteration a
  mesh smoothing is performed.}

  \item{smoothit}{integer: determine smoothing iterations.}

  \item{smoothtype}{character: select smoothing algorithm -
  see vcgSmooth for further details.}

  \item{tol}{numeric: convergence threshold of MSE between
  vertices of two successive iterations.}

  \item{useiter}{logical: if TRUE, each iteration uses the
  updated reference mesh, if false. The original mesh will
  be deformed based on the updated correspondences.}

  \item{minclost}{minimum amount of correspondence points.
  If less correspondences are found, dist will be increased
  by "distinc" (see below).}

  \item{distinc}{increment of dist, in case minclost is not
  reached.}

  \item{scale}{logical: if TRUE, initial landmark based
  rigid registration includes scaling.}

  \item{icp}{vector of length 4. Passing parameters to
  \code{\link{icp}}, which is performed after intial
  landmark based registration. The parameters are
  icp[1]=iterations; icp[2]=rhotol; icp[3]=uprange, and
  icp[4]=scale. If icp=NULL, no ICP-matching is performed.
  E.g. icp=c(3,pi/2,0.6,TRUE) will result in 3 icp
  iterations, condidering the closest 60\% of
  correspondences with normal deviation of pi/2 and include
  scaling.}
}
\value{
  \item{mesh}{registered mesh} \item{meshrot }{mesh1,
  rotated onto mesh2} \item{lm1rot }{lm1, rotated onto lm2}
  \item{lmtmp1 }{correspondences on updated reference mesh
  of last iteration} \item{lmtmp2 }{correspondences on
  updated target mesh of last iteration}
}
\description{
  Perform registration of two triangular meshes, minimizing
  per-face distortions.
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

