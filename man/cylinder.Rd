\name{cylinder}
\alias{cylinder}
\title{create a mesh representation of a cylinder}
\usage{
cylinder(x, dirs, length, radius = 1, fine = 20, adNormals = FALSE)
}
\arguments{
  \item{x}{coordinates of center of the cylinders base -
  either a 3D-vector of a matrix of dim k x 3}

  \item{dirs}{directions of the cylinder to point - either
  a 3D-vector of a matrix of dim k x 3}

  \item{length}{length of the cylinder}

  \item{radius}{radius of the cylinder's base}

  \item{fine}{amount of vertices to create the base circle}

  \item{adNormals}{logical: if TRUE, normal vectors for
  each vertex are calculated.}
}
\value{
a triangular mesh of class 'mesh3d'
}
\description{
create a mesh representation of a cylinder
}
\examples{
require(rgl)
wire3d(cylinder(1:3,1:3,4,adNormals = TRUE))
}

