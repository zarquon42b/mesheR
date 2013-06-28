\name{normcheck}
\alias{normcheck}
\title{compare normal directions between two states of a mesh}
\usage{
  normcheck(mesh1, mesh2, circle = TRUE)
}
\arguments{
  \item{mesh1}{triangular mesh}

  \item{mesh2}{triangular mesh}

  \item{circle}{logical: which method to use calculating
  the angel}
}
\value{
  numeric vector containing angles between corresponding
  normals
}
\description{
  compare normal directions between two states of a mesh
}
