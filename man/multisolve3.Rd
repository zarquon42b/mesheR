\name{multisolve3}
\alias{multisolve3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fortran routine for fast inversion of multiple 3x3 matrices
}
\description{Fortran routine for fast inversion of multiple 3x3 matrices

 }
\usage{
multisolve3(A, trans = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{m x 3 matrix, where m is a multiple of 3
%%     ~~Describe \code{A} here~~
}
  \item{trans}{calculate inverse of transposed subblocks
%%     ~~Describe \code{trans} here~~
}
}
\details{Be m=k*3. The matrix will then be interpreted as k 3x3
  matrices, which are subsequently inverted
  
}
\value{Returns matrix containing inverted submatrices.

}

}
\author{Stefan Schlager
}


\examples{
A <- matrix(rnorm(18),6,3)
B <- multisolve3(A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
