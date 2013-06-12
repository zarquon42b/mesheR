\name{AmbergRegister}
\alias{AmbergRegister}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Register two triangular meshes based on smooth deformation.
}
\description{Perform registration of two triangular meshes, minimizing
  per-face distortions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AmbergRegister(mesh1, mesh2, lm1 = NULL, lm2 = NULL, k = 1, lambda = 1, iterations = 15, rho = pi/2, dist = 2, border = FALSE, smooth = TRUE, tol = 1e-04, useiter = TRUE, minclost = 50, distinc = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh1}{reference mesh: triangular mesh of class "mesh3d". No loose vertices, edges and
    degenerated faces are allowed.
%%     ~~Describe \code{mesh1} here~~
}
  \item{mesh2}{target mesh: triangular mesh of class "mesh3d". 
%%     ~~Describe \code{mesh2} here~~
}
  \item{lm1}{m x 3 matrix containing correspondences on "mesh1"
%%     ~~Describe \code{lm1} here~~
}
  \item{lm2}{m x 3 matrix containing target correspondences on "mesh2"
%%     ~~Describe \code{lm2} here~~
}
  \item{k}{integer: parameter regularizing face normal distortion. Can
 be vector of length(iterations) or single value.
%%     ~~Describe \code{k} here~~
}
  \item{lambda}{numeric: parameter regularizing faces's distortion. Can
 be vector of length(iterations) or single value.
%%     ~~Describe \code{lambda} here~~
}
  \item{iterations}{integer: number of iterations to run.
%%     ~~Describe \code{iterations} here~~
}
  \item{rho}{numeric: 0 < rho < 2*pi tolerance of normal deviation between reference vertices
    and corresponding closest points on target suface.
%%     ~~Describe \code{rho} here~~
}
  \item{dist}{numeric: tolerance of maximal distance between reference vertices
    and corresponding closest points on target suface.
%%     ~~Describe \code{dist} here~~
}
  \item{border}{logical: if FALSE, hits on border faces are ignored
    (reduces distortion)
%%     ~~Describe \code{border} here~~
}
  \item{smooth}{logical: if TRUE after each iteration a Taubin smooth is
    performed (1 iteration).
%%     ~~Describe \code{smooth} here~~
}
  \item{tol}{numeric: convergence threshold of MSE between vertices of
    two successive iterations.
%%     ~~Describe \code{tol} here~~
}
  \item{useiter}{logical: if TRUE, each iteration uses the updated
    reference mesh, if false. The original mesh will be deformed based
    on the updated correspondences.
%%     ~~Describe \code{useiter} here~~
}
  \item{minclost}{minimum amount of correspondence points. If less
    correspondences are found, dist will be increased by "distinc" (see below).
%%     ~~Describe \code{minclost} here~~
}
  \item{distinc}{increment of dist, in case minclost is not reached.
%%     ~~Describe \code{distinc} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 return(list(mesh = tmp$mesh, meshrot = meshorig, lm1rot = lm1, 
        lmtmp1 = lmtmp1, lmtmp2 = lmtmp2))%%  ~Describe the value returned

  \item{mesh }{registered mesh}
  \item{meshrot }{mesh1, rotated onto mesh2}
  \item{lm1rot }{lm1, rotated onto lm2}
  \item{lmtmp1 }{correspondences on updated reference mesh of last iteration}
  \item{lmtmp2 }{correspondences on updated target mesh of last iteration}

%% ...
}
\references{
Amberg, B. 2011. Editing faces in videos, University of Basel.
}
\author{
 Stefan Schlager
}


\seealso{
\code{\link{gaussMatch}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
