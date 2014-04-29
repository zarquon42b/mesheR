\name{pPCAcond}
\alias{pPCAcond}
\alias{restrictMissing}
\alias{setMod}
\alias{setMod.pPCAcond}
\title{calculate or modify a conditional PCA based on 3D-coordinates}
\usage{
pPCAcond(array, missingIndex, sigma = NULL, exVar = 1)

setMod(procMod, sigma, exVar)

\method{setMod}{pPCAcond}(procMod, sigma = NULL, exVar = 1)

restrictMissing(x, model)
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

  \item{cores}{integer how many cores to use.}
}
\value{
return a probabilistic PCA model of class "procMod"
}
\description{
calculate or modify a conditional PCA based on
3D-coordinates
}
\examples{
require(Morpho)
data(boneData)
model <- pPCAcond(boneLM[,,-1],missingIndex=3:4)
## change parameters without recomputing Procrustes fit
model1 <- setMod(model, sigma=1, exVar=0.8)
}

