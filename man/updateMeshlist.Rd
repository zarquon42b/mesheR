% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updateMeshlist.r
\name{updateMeshlist}
\alias{updateMeshlist}
\title{replace vertices in a meshlist with new values from an array}
\usage{
updateMeshlist(x, array)
}
\arguments{
\item{x}{list containing triangular meshes of class mesh3d}

\item{array}{k x 3 x n array containing replacement vertices for x; the i-th matrix will be inserted in the i-th mesh}
}
\value{
list with updated meshes
}
\description{
replace vertices in a meshlist with new values from an array
}
