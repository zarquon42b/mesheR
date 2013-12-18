\name{placeDowels}
\alias{placeDowels}
\title{Place Dowels representing tissue depth on a skull}
\usage{
placeDowels(lm, mesh, ldowel, render = TRUE, col = 1, radius = 1,
  meshcol = 3, fine = 50, smooth = TRUE, dowelcol = NULL)
}
\arguments{
  \item{lm}{k x 3landmarks representing dowel locations}

  \item{mesh}{triangular mesh (object of class 'mesh3d')}

  \item{ldowel}{vector of length k giving a length for each
  dowel}

  \item{render}{logical: if TRUE, the result will be
  rendered in an rgl window}

  \item{col}{color of the dowels}

  \item{radius}{numeric: diameter of dowels}

  \item{meshcol}{mesh color}

  \item{fine}{integer: amount of vertices to generate the
  cylinders representing tissue thickness.}

  \item{smooth}{logical: use smoothed normals for dowel
  orientation}

  \item{dowelcol}{specify dowelcolor.}
}
\value{
an invisible list with \item{dowels }{a list containing the
meshes for all dowels} \item{endpoints }{a matrix
containing the coordinates of the dowels endpoints}
}
\description{
Place Dowels representing tissue depth on a skull by using
predefined tissue-values
}
\details{
For orientation of the dowels, the (angle weighted) normal
vectors of the surface is used.
}
\examples{
require(Rvcg)
data(humface)
lms <- matrix(c(17.6061 , 9.2072 , -6.9917 , 44.7959 , 37.1135 , 76.0469 , 3.7734 , 12.8234 , 81.3138 , 63.4207 , 47.5009 , 44.6468), 4, 3)
#place markers of 5mm length
placeDowels(lms, humface, rep(5, 4),output = FALSE)
}
\seealso{
\code{\link{cylinder}}
}

