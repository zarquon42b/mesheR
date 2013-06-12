\name{AmbergDeformSpam}
\alias{AmbergDeformSpam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Deform triangular mesh based on correspondences
}
\description{Perform smooth deformation of a triangular mesh, minimizing
  per-face distortions.
}
\usage{
AmbergDeformSpam(mesh, lm1, lm2, k0 = 1, lambda = 1, S = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{
    triangular mesh of class "mesh3d". No loose vertices, edges and
    degenerated faces are allowed.
  }
  \item{lm1}{m x 3 matrix containing correspondences on "mesh"
  }
  \item{lm2}{m x 3 matrix containing target correspondences
    %%     ~~Describe \code{lm2} here~~
  }
  \item{k0}{integer: parameter regularizing face normal distortion.
  }
  \item{lambda}{numeric: parameter regularizing faces's distortion.
      %%     ~~Describe \code{lambda} here~~
    }
    \item{S}{optional: object from function createS from previous calculation.
      %%     ~~Describe \code{S} here~~
    }
  }
  \details{Perform smooth deformation of a triangular mesh, minimizing
    per-face distortions.No loose vertices, edges and
    degenerated faces are allowed, as they lead to singular equation system.
  }
  \value{
    \item{mesh }{deformed mesh}
    \item{Jn }{Jacobi submatrix Jn}
    \item{Jc }{Jacobi submatrix Jc}
    \item{J }{Jacobian matrix}
    \item{H }{Hessian of J, class "spam"}
    \item{Hchol }{Cholesky decomposition of H; class"spam"}

    %% ...
  }
  \references{Amberg, B. 2011. Editing faces in videos, University of Basel.
  }
  \author{
    Stefan Schlager}


  \seealso{
    \code{\link{gaussMatch}}
  }

  % Add one or more standard keywords, see file 'KEYWORDS' in the
  % R documentation directory.
  \keyword{ ~kwd1 }
  \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
