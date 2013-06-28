\name{icp}
\alias{icp}
\title{Iterative closest point matching between two triangular meshes.}
\usage{
  icp(mesh1, mesh2, iterations = 3, scale = T, lm1 = NULL,
    lm2 = NULL, uprange = 0.9, rhotol = pi, k = 50,
    reflection = FALSE, pro = c("morpho", "vcg"))
}
\arguments{
  \item{mesh1}{object of class "mesh3d".}

  \item{mesh2}{object of class "mesh3d".}

  \item{iterations}{integer}

  \item{scale}{lgoical: scale reference mesh to optimize
  fit.}

  \item{lm1}{optional: kx3 matrix containing reference
  points on mesh1.}

  \item{lm2}{optional: kx3 matrix containing reference
  points on mesh2.}

  \item{uprange}{quantile of distances between vertices of
  mesh1 and closest points on mesh2. All hit points on
  mesh2 farther away than the specified quantile are not
  used for matching.}

  \item{rhotol}{maximum allowed angle of which normals are
  allowed to differ between reference points and hit target
  points. Hit points above this threshold will be
  neglected.}

  \item{k}{integer: number of closest triangles on target
  that will be searched for hits. Only used when
  pro="morpho".}

  \item{reflection}{logical: allow reflection.}

  \item{pro}{character: algorithm for closest point search
  "morpho" calls closemeshKD from the package Morpho, while
  vcg calls vcgClost from mesheR. If the region of the
  targetmesh is much smaller than the region of the
  reference, "vcg" can be really slow. Otherwise very fast.
  "morpho" is the stable but somewhat slower algorithm.}
}
\value{
  Returns the rotated mesh1.
}
\description{
  performs rigid body transformations (and scaling if
  requested) to map a reference mesh onto a target mesh
}
\details{
  the registration is done by minimising squared distances
  between reference vertices and closest points on the
  target surface (a.k.a. Procrustes registration)
}
\examples{
data(nose)
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
rotnose <- icp(warpnose.long,shortnose.mesh,lm1=longnose.lm,lm2=shortnose.lm,rhotol=0.7,uprange=0.9)
shade3d(rotnose,col=2,alpha=0.7)
shade3d(shortnose.mesh,col=3,alpha=0.7)
}
\author{
  Stefan Schlager
}
\references{
  Zhang Z. 1994. Iterative point matching for registration
  of free-form curves and surfaces International Journal of
  Computer Vision 13:119-152.
}
\seealso{
  \code{\link{rotmesh.onto}}, \code{\link{rotonto}}
}
\keyword{~kwd1}
\keyword{~kwd2}

