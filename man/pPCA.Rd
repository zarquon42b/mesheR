\encoding{utf8}
\name{pPCA}
\alias{pPCA}
\alias{pPCAcond}
\alias{predictPCAcond}
\alias{predictPCAcond.matrix}
\alias{predictPCAcond.mesh3d}
\alias{predictpPCA}
\alias{predictpPCA.matrix}
\alias{predictpPCA.mesh3d}
\alias{predictpPCA.numeric}
\alias{setMod}
\alias{setMod.pPCA}
\alias{setMod.pPCAcond}
\title{calculate or modify a probablistic PCA based on 3D-coordinates}
\usage{
pPCA(array, sigma = NULL, exVar = 1, scale = TRUE, refmesh = NULL)

pPCAcond(array, missingIndex, sigma = NULL, exVar = 1, refmesh = NULL,
  scale = TRUE)

setMod(procMod, sigma, exVar)

\method{setMod}{pPCA}(procMod, sigma = NULL, exVar = 1)

\method{setMod}{pPCAcond}(procMod, sigma = NULL, exVar = 1)

predictPCAcond(x, model, refmesh, sdmax, pPCA = FALSE, ...)

\method{predictPCAcond}{matrix}(x, model, refmesh = FALSE, sdmax,
  origSpace = TRUE, pPCA = FALSE, ...)

\method{predictPCAcond}{mesh3d}(x, model, refmesh = FALSE, sdmax,
  origSpace = TRUE, pPCA = FALSE, ...)

predictpPCA(x, model, refmesh = FALSE, ...)

\method{predictpPCA}{matrix}(x, model, refmesh = FALSE, sdmax = 2,
  origSpace = TRUE, ...)

\method{predictpPCA}{mesh3d}(x, model, refmesh = FALSE, sdmax = 2,
  origSpace = TRUE, ...)

\method{predictpPCA}{numeric}(x, model, refmesh = FALSE, ...)
}
\arguments{
  \item{array}{array of dimensions k x 3 x n, where
  k=number of coordinates and n=sample size.}

  \item{missingIndex}{integer vector: specifies which
  points are missing in the conditional model}

  \item{procMod}{object of class "procMod" as returned by
  pPCAcond or setMod}

  \item{sigma}{estimate error variance (sensible is a value
  estimating coordinate error in terms of observer error)}

  \item{exVar}{numeric value with \code{0 < exVar <= 1}
  specifying the PCs to be included by their cumulative
  explained Variance}

  \item{refmesh}{a triangular mesh, where the vertices
  correspond to the coordinates in \code{array}}

  \item{scale}{logical: allow scaling in Procrustes
  fitting}

  \item{x}{vector of deviation in standard deviations,
  coordinate matrix or triangular mesh of class "mesh3d" to
  be predicted.}

  \item{sdmax}{maximum allowed standard deviation (per
  Principal axis) within the model space. Defines the
  probabilistic boundaries.}

  \item{origSpace}{logical: rotate the estimation back into
  the original coordinate system.}

  \item{pPCA}{logical: if TRUE, a constrained pPCA model is
  returned.}

  \item{model}{probabilistic model of class "pPCA" or
  "pPCAcond"}
}
\value{
\code{pPCA} and \code{pPCAcond} return a probabilistic PCA
model of class "pPCA" or "pPCAcond" respectively.
\code{predictPCA} and \code{predictPCAcond} select the most
probable shape within a given model (within defined
boundaries), \code{setMod} is used to modify existing
models without recomputing Procrustes registration and
subsequent PCA.
}
\description{
calculate or modify a probablistic PCA based on
3D-coordinates

#'

#'
}
\examples{
require(Morpho)
data(boneData)
model <- pPCAcond(boneLM[,,-1],missingIndex=3:4)
## change parameters without recomputing Procrustes fit
model1 <- setMod(model, sigma=1, exVar=0.8)
}
\references{
\enc{Lüthi}{Luethi} M, Albrecht T, Vetter T. 2009.
Probabilistic modeling and visualization of the flexibility
in morphable models. In: Mathematics of Surfaces XIII.
Springer. p 251-264
}

