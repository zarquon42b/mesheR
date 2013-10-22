\name{AmbergRegister}
\alias{AmbergRegister}
\title{Register two triangular meshes based on smooth deformation.}
\usage{
  AmbergRegister(mesh1, mesh2, lm1 = NULL, lm2 = NULL,
    k = 1, lambda = 1, iterations = 15, rho = pi/2,
    dist = 2, border = FALSE, smooth = TRUE, smoothit = 1,
    smoothtype = "t", tol = 1e-04, useiter = TRUE,
    minclost = 50, distinc = 1, scale = TRUE,
    reflection = FALSE, icp = NULL, nn = 20, cores = 1,
    silent = FALSE)
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

  \item{reflection}{logical: if TRUE, initial landmark
  based rigid registration allows reflections.}

  \item{icp}{vector of length 4. Passing parameters to
  \code{\link{icp}}, which is performed after intial
  landmark based registration. The parameters are
  icp[1]=iterations; icp[2]=rhotol; icp[3]=uprange, and
  icp[4]=scale. If icp=NULL, no ICP-matching is performed.
  E.g. icp=c(3,pi/2,0.6,TRUE) will result in 3 icp
  iterations, condidering the closest 60\% of
  correspondences with normal deviation of pi/2 and include
  scaling.}

  \item{nn}{integer: closest barycenters. During search for
  closest points on target, the closest \code{nn} faces are
  probed. The larger \code{nn} is , the more accurate the
  closest point search but also the more time consuming.}

  \item{cores}{integer: how many cores to use for closest
  point search}

  \item{silent}{logical: no verbosity}
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
\details{
  This function runs an elastic-ICP surface matching
  algorithm, that minimizes the original meshes internal
  structure by solving a sparse equation system. The user
  can control 2 parameters of mesh stiffness: \code{lambda}
  and \code{k}. \code{lambda} controls the impact of the
  control points (closest points) as it is a weight applied
  to the equation system. The value of \code{lambda} should
  be carefully selected depending on the object overall
  size: i.e. to match two tiny meshes one will need a
  higher value than a for a larger object (example: I found
  values between 0 and 1 suitable for human faces and
  values between 10 and 100 suitable for mice teeth).
  \code{k} controls the normal slackness, i.e. the
  deviation of normal direction. The larger, \code{k}, the
  more elastic the deformation will be. \code{lambda} and
  \code{k} can be specified as vectors of length
  \code{iterations}, to assign a specific value for each
  iteration.
}
\examples{
require(Morpho)
require(rgl)
require(Rvcg)
data(humface)
data(dummyhead)
## set parameters making each iteration more elastic
# only 10 iterations to keep example calculation time reasonable.
params <- list(iterations=10)
params <- append(params, list(
   # first \\code{lambda} is set relatively high because first matching uses landmarks
   # then let it increase from 0.2 to 0.6
   lambda=c(0.7,seq(from = 0.2,to=0.6,length.out = params$iterations-1)),
   # treat \\code{k} similar as \\code{lambda}
   k=c(10,seq(from = 1,to=params$iterations-1,by=1)),
   useiter=FALSE # iteratively deform dummyhead onto humface
   ))
map <- AmbergRegister(dummyhead.mesh, humface, lm1=dummyhead.lm,
                 lm2=humface.lm, iterations=params$iterations,
                 k=params$k, lambda=params$lambda, useiter=params$useiter)
# compare matched and original face:
meshDist(map$mesh, humface ,from=-3,to=3,tol=0.5)
# render original mesh as wireframe
wire3d(humface)
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

