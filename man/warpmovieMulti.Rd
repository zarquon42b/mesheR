\name{warpmovieMulti}
\alias{warpmovieMulti}
\title{make a warpmovie with multiple registered meshes}
\usage{
  warpmovieMulti(..., n, col = "green", folder = NULL,
    movie = "warpmovie", add = FALSE, close = TRUE,
    countbegin = 0, ask = TRUE, whichcolor = 1,
    align = TRUE, scale = FALSE)
}
\arguments{
  \item{\dots}{registered meshes}

  \item{n}{amount of intermediate images between two
  meshes}

  \item{col}{mesh color}

  \item{folder}{character: output folder for created images
  (optional)}

  \item{movie}{character: name of the output files}

  \item{add}{logical: if TRUE, the movie will be added to
  the focussed rgl-windows.}

  \item{close}{logical: if TRUE, the rgl window will be
  closed when finished. width and 200 the height of the
  image.}

  \item{whichcolor}{select which mesh's vertexcolors to
  use.}

  \item{align}{logical: align all meshes to the first one}

  \item{scale}{if \code{align=TRUE} this controls if
  scaling is allowed.}

  \item{countbegin}{integer: number to start image
  sequence.}

  \item{ask}{logical: if TRUE, the viewpoint can be
  selected manually.}
}
\description{
  make a warpmovie with multiple registered meshes - this
  is a pimped version of warpmovie3d from package 'Morpho'
}

