% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/displacementField.r
\name{applyDisplacementField}
\alias{applyDisplacementField}
\title{apply a discrete displacement field to a set of points/mesh in its domain}
\usage{
applyDisplacementField(
  dispfield,
  points,
  k = 10,
  sigma = 20,
  type = c("Gauss", "Laplace", "Exponential", "TPS"),
  gamma = 1,
  lambda = 1e-08,
  threads = 1
)
}
\arguments{
\item{dispfield}{displacement field of class "DisplacementField", e.g. created using \code{\link{createDisplacementField}} or}

\item{points}{matrix or mesh3d at which to evaluate the interpolated displacement field}

\item{k}{integer: number of k closest points to evaluate.}

\item{sigma}{kernel bandwidth used for smoothing. For all kernels except B-spline, sigma controls the importance of the neighbourhood by defining the bandwidth of the smoothing kernel. For B-spline it defines the support (the higher, the "wobblier" the deformation field can become.}

\item{type}{kernel function for smoothing are "Gauss","Laplace", "Exponential" and "Bspline" (or any abbreviation thereof).}

\item{gamma}{dampening factor (displacement vectors will be divided by \code{gamma}}

\item{lambda}{smoothing factor for TP}

\item{threads}{integer: number of threads to use for computing the interpolation.}
}
\value{
returns the displaced version of points
}
\description{
apply a discrete displacement field to a set of points/mesh in its domain by applying the gaussian smoothed interpolation based of k closest neighbours
}
\note{
if points is identical to the domain of the displacement field, no interpolation will be performed.
}
