\name{gaussMatch}
\alias{gaussMatch}
\title{map two surface meshes using smoothed displacement fields}
\usage{
  gaussMatch(mesh1, mesh2, iterations = 10, smooth = NULL,
    smoothit = 10,
    smoothtype = c("taubin", "laplace", "HClaplace"),
    sigma = 20, gamma = 2, f = 1.2, oneway = F, lm1 = NULL,
    lm2 = NULL, icp = FALSE, icpiter = 3, uprange = 0.95,
    rhotol = 1, nh = NULL, toldist = 0, patch = NULL,
    repro = FALSE, cores = detectCores(),
    pro = c("morpho", "vcg"), k0 = 50, prometh = 1,
    angtol = NULL, border = FALSE, horiz.disp = NULL,
    Amberg = NULL, silent = FALSE, ...)
}
\arguments{
  \item{mesh1}{An object of class mesh3d used as atlas mesh
  to be deformed onto the target mesh. Mesh resolution
  should be ~1.5.}

  \item{mesh2}{An object of class mesh3d used as target
  mesh. Mesh resolution should be ~1.5.}

  \item{iterations}{Iterations of displacement. Default is
  10.}

  \item{smooth}{Integer: smoothing factor. Default is NULL,
  no smoothing.}

  \item{smoothit}{integer: smoothing steps.}

  \item{smoothtype}{Type of smoothing: Taubin, Laplacian,
  or HClaplacian. For details see \code{\link{vcgSmooth}}}

  \item{sigma}{starting parameter for smoothed displacement
  (see Moshfeghi 1994). Sigma controls the importance of
  the neighbourhood by defining the standard-deviation for
  the gaussian smoothing}

  \item{gamma}{stiffness factor controlling displacement
  strength. The smoothed displacement vector for each
  vertex is divided by \code{gamma}. The larger
  \code{gamma}, the slower the approximation.}

  \item{f}{parameter controlling iterative decrease of
  \code{sigma} making the displacement locally more elastic
  with each iteration. (Moshfeghi 1994). Starting with
  \code{sigma}, this parameter for the k-th iteration is
  \code{sigma *f ^(-k)}}

  \item{oneway}{logical: only displace towards the target
  without taking into account the displacement from the
  target.}

  \item{lm1}{A k x 3 matrix containing landmarks
  corrresponding to mesh1 for initial rotation of mesh1
  onto mesh2.}

  \item{lm2}{A k x 3 matrix containing landmarks
  corrresponding to mesh2 for initial rotation of mesh1
  onto mesh2.}

  \item{icp}{Logical: if TRUE, iterative closest point
  procedure will be executed.}

  \item{icpiter}{Integer: Number of iterations of icp.}

  \item{uprange}{argument passed to icp (see
  \code{\link{icp}})}

  \item{rhotol}{Numeric: argument passed to
  \code{\link{icp}}.Exclude target points with deviation of
  normals larger than than rhotol.}

  \item{nh}{Integer: neighbourhood (number vertices) for
  controlling displacement smoothing, default is 150/mesh
  resolution.}

  \item{toldist}{Integer: Exclude everything from the whole
  procedure with a greater distance from initial point than
  toldist. 0 disables this feature.}

  \item{patch}{A m x 3 matrix containing the atlas landmark
  configuration on mesh1 which is not present in mesh2 and
  will automatically placed on the latter.}

  \item{repro}{Logical: If TRUE, a reprojection of patch
  onto the iteratively estimated target surface will be
  performed after each iteration.}

  \item{cores}{integer: number of CPU cores to be used for
  multithreaded subroutines.}

  \item{pro}{which projection method to use: "m"=
  \code{\link{closemeshKD}} from Morpho; "v"=
  \code{\link{vcgClost}} from package Rvcg}

  \item{k0}{Integer: argument passed to closemeshKD (will
  be argument "k" in \code{\link{closemeshKD}} .}

  \item{prometh}{argument passed to closemeshKD.  Integer:
  0 or 1. If prometh=0, take closest point for
  displacement. If prometh=1, do not just take the closest
  point, but for two absolut distances which are the same,
  take the point which is orthogonal to the closest face
  see Moshfeghi 1994).}

  \item{angtol}{numeric: If the angle between hit points'
  normals and the starting points' normals exceeds this
  threshold the displacement vector will be discarded.}

  \item{border}{Logical: if TRUE, displacement vectors
  hitting mesh borders are discarded.}

  \item{horiz.disp}{numeric: If the angle between hit
  points' normals (independent of its orientation) and the
  distance vector between hit point and starting points
  exceeds this threshold, the displacement vector will be
  discarded. Reduces distortion especially at mesh
  borders.}

  \item{Amberg}{vector containing 2 arguments invoking a
  smooth Deformation using \code{\link{AmbergDeformSpam}}.
  Layout: Amgerg=c(lambda, k0)}

  \item{silent}{logical suppress messages}

  \item{\dots}{Further arguments passed to \code{nn2}.}
}
\value{
  If a patch is specified: \item{mesh}{matched mesh}
  \item{patch}{displaced patch as specified in input.} else
  a mesh of class "mesh3d" is returned.
}
\description{
  Map a reference mesh onto a target surface using
  displacement fields.
}
\details{
  This function implements the mesh matching method
  suggested by Moshfeghi et al. and Bryan et al..
  Additional mechanisms for controlling and restricting the
  displacement smoothing are implemented
}
\examples{
require(Morpho)
data(nose)##load data
##warp a mesh onto another landmark configuration:
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
### result won't be too good as the surfaces do stronly differ.
match <- gaussMatch(shortnose.mesh,warpnose.long,gamma=4,iterations=3,smooth=1,smoothtype="h",smoothit=10,nh=50,angtol=pi/2)
}
\author{
  Stefan Schlager
}
\references{
  Bryan, R., Mohan, P. S., Hopkins, A., Galloway, F.,
  Taylor, M., and Nair, P. B. 2010. Statistical modelling
  of the whole human femur incorporating geometric and
  material properties. Medical Engineering &amp; Physics,
  32(1):57-65.

  Moshfeghi, M., Ranganath, S., and Nawyn, K. 1994.
  Three-dimensional elastic matching of volumes. IEEE
  Transactions on Image Processing: A Publication of the
  IEEE Signal Processing Society, 3(2):128-138.
}
\seealso{
  \code{\link{meshres}}, \code{\link{vcgClost}},
  \code{\link{vcgBorder}}, \code{\link{icp}},
  \code{\link{vcgSmooth}}
}
\keyword{~kwd1}
\keyword{~kwd2}

