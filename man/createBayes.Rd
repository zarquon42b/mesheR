\name{createBayes}
\alias{createBayes}
\title{create an object of class BayesDeform to be used in gaussMatch}
\usage{
createBayes(model, sd = 3, nPC = NULL, scale = FALSE, maxVar = 95)
}
\arguments{
  \item{model}{list containing mshape, PCs, sd of a model
  distribution (can be result of \code{procSym}}

  \item{sd}{numeric or vector of numerics (length must be
  of value =< \code{iterations} as to be specified in
  gaussMatch or AmbergRegister. If \code{length(sd) =n <
  iterations}, then the first n iterations will be
  restricted to the according standard deviations}

  \item{nPC}{number of PCs used in restricting (overrides
  maxVar)}

  \item{scale}{logical: if FALSE, configuration will be
  restricted within a probablistic hypercube, otherwise
  Chi-Square distribution will be used.}

  \item{maxVar}{percentage to be explained by the PCs used
  for the restriction. E.g. if you want to include PCs
  explaining 95\% of the variance, set maxVar=95.}
}
\value{
returns a list of class "BayesDeform"
}
\description{
create an object of class BayesDeform to be used in
gaussMatch
}

